package org.broadinstitute.hellbender.tools.spark.sv.prototype;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariationSparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.DiscoverVariantsFromContigAlignmentsSAMSpark;
import org.broadinstitute.hellbender.tools.spark.sv.utils.RDDUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;

/**
 * Tool to select alignments of long reads possibly pointing to pathogen integration site.
 */
@CommandLineProgramProperties(summary="Tool to select alignments of long reads possibly pointing to pathogen integration site.",
        oneLineSummary="Tool to select alignments of long reads possibly pointing to pathogen integration site.",
        omitFromCommandLine = true,
        usageExample = "gatk-launch ExtractSuspectedPathogenAlignmentsSpark \\" +
                "-I <PATH_TO_ASSEMBLY_ALN.sam> \\" +
                "-O <OUTPUT.sam> \\" +
                "-R <PATH_TO_REF> \\" +
                "--maxGapSize 150",
        programGroup = StructuralVariationSparkProgramGroup.class)
@BetaFeature
public class ExtractSuspectedPathogenAlignmentsSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;
    private final Logger localLogger = LogManager.getLogger(ExtractSuspectedPathogenAlignmentsSpark.class);

    @Argument(doc = "sam file for aligned contigs", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String outputAssemblyAlignments;

    @Argument(doc = "size of the maximum gap (uncovered by alignments) a read must have for it to be included in output",
            shortName = "mgs", fullName = "maxGapSize")
    private int maxGapSize;

    @Override
    public boolean requiresReads() {
        return true;
    }

    @Override
    public boolean requiresReference() {
        return true;
    }


    @Override
    protected void runTool(final JavaSparkContext ctx) {

        final JavaRDD<GATKRead> rawAlignments = getUnfilteredReads();
        final Tuple2<JavaRDD<GATKRead>, JavaRDD<GATKRead>> unmappedAndMapped = RDDUtils.split(rawAlignments, GATKRead::isUnmapped, true);
        final JavaRDD<GATKRead> unmappedReads = unmappedAndMapped._1;
        localLogger.info("Found " + unmappedReads.count() + " unmapped contigs.");

        final HashSet<String> selectedMappedContigNames =
                new HashSet<>(selectMappedContigs(unmappedAndMapped._2, maxGapSize, getHeaderForReads(), localLogger)
                              .map(tig -> tig.contigName).distinct().collect());
        writeReads(ctx, outputAssemblyAlignments, rawAlignments.filter(r -> selectedMappedContigNames.contains(r.getName())).union(unmappedReads));
    }

    static JavaRDD<AlignedContig> selectMappedContigs(final JavaRDD<GATKRead> mappedReads, final int coverageThresholdInclusive,
                                                      final SAMFileHeader header, final Logger toolLogger) {

        final JavaRDD<AlignedContig> filter = new DiscoverVariantsFromContigAlignmentsSAMSpark
                .SAMFormattedContigAlignmentParser(mappedReads, header, true, toolLogger)
                .getAlignedContigs()
                .filter(ctg -> keepContigForPathSeqUse(ctg, coverageThresholdInclusive));
        toolLogger.info("Found " + filter.count() + " mapped contigs suggesting pathogen injection sites");
        return filter;
    }

    /**
     * Currently keep contigs that has
     *   1) alignments covering less than or equal to only half of the contig sequence,
     *   3) a maximum gap (uncovered by any alignment) size that is large than or equal to {@code maxGapSzThresholdInclusive}
     */
    @VisibleForTesting
    static boolean keepContigForPathSeqUse(final AlignedContig contig,
                                           final int maxGapSzThresholdInclusive) {

        final Tuple2<Integer, Integer> coverageAndMaxGapSz = alignmentsCoverage(contig);
        final int alignmentCoverage = coverageAndMaxGapSz._1,
                  maxGapSz = coverageAndMaxGapSz._2;

        return contig.contigSequence.length >= 2 * alignmentCoverage || maxGapSz >= maxGapSzThresholdInclusive;
    }

    /**
     * @return  a paired value of (how long of the provided contig is covered with alignments,
     *                             largest uncovered gap size)
     */
    @VisibleForTesting
    static Tuple2<Integer, Integer> alignmentsCoverage(final AlignedContig contig) {
        if (contig.alignmentIntervals.isEmpty()) return DEFAULT_EMPTY_COV;

        final List<AlignmentInterval> alignmentIntervals = contig.alignmentIntervals;
        final List<SVInterval> maximallyExtendedCovers = new ArrayList<>(alignmentIntervals.size());

        final Iterator<AlignmentInterval> it = alignmentIntervals.iterator();
        AlignmentInterval current = it.next(); // at least one exists
        final int contigLength = contig.contigSequence.length;
        int maxUncoveredGapSize = Math.max(current.startInAssembledContig - 1,
                                           contigLength - alignmentIntervals.get(alignmentIntervals.size()-1).endInAssembledContig); // head and tail
        SVInterval currentSVI = new SVInterval(1, current.startInAssembledContig, current.endInAssembledContig+1); // +1 for [a,b) semi-closed
        while(it.hasNext()) {
            final AlignmentInterval next = it.next();
            final SVInterval nextSVI = new SVInterval(1, next.startInAssembledContig, next.endInAssembledContig+1);

            if (nextSVI.overlaps(currentSVI)) {
                currentSVI = new SVInterval(1, Math.min(currentSVI.getStart(), nextSVI.getStart()),
                                                      Math.max(currentSVI.getEnd(), nextSVI.getEnd()));
            } else {
                maximallyExtendedCovers.add(currentSVI);
                currentSVI = nextSVI;
                maxUncoveredGapSize = Math.max(maxUncoveredGapSize,
                                               next.startInAssembledContig - current.endInAssembledContig - 1);
            }
        }
        maximallyExtendedCovers.add(currentSVI);

        return new Tuple2<>(maximallyExtendedCovers.stream().mapToInt(SVInterval::getLength).sum(), maxUncoveredGapSize);
    }

    private static final Tuple2<Integer, Integer> DEFAULT_EMPTY_COV = new Tuple2<>(0,0);
}