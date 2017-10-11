package org.broadinstitute.hellbender.tools.spark.sv.prototype;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.collections.iterators.EmptyListIterator;
import org.apache.commons.collections.iterators.SingletonListIterator;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariationSparkProgramGroup;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.*;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.prototype.SvDiscoverFromLocalAssemblyContigAlignmentsSpark;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.AlignedAssemblyOrExcuse;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.FindBreakpointEvidenceSpark;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVVCFWriter;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import scala.Tuple2;

import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;

import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection.DEFAULT_MIN_ALIGNMENT_LENGTH;

/**
 * Tool to run the sv pipeline up for possible pathogen integration site detection and assembled contigs with low alignment coverage.
 */
@CommandLineProgramProperties(summary="Tool to run the sv pipeline up for possible pathogen integration site detection and assembled contigs with low alignment coverage.",
        oneLineSummary="Tool to run the sv pipeline up for possible pathogen integration site detection and assembled contigs with low alignment coverage.",
        omitFromCommandLine = true,
        usageExample = "gatk-launch DetectPossiblePathogenIntegraionBkptPipelineSpark \\" +
                "-I <PATH_TO_ASSEMBLY_ALN.sam> \\" +
                "-O <OUTPUT.sam> \\" +
                "-R <PATH_TO_REF> \\" +
                "--maxGapSize 150",
        programGroup = StructuralVariationSparkProgramGroup.class)
@BetaFeature
public class DetectPossiblePathogenIntegraionBkptPipelineSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;
    private final Logger localLogger = LogManager.getLogger(DetectPossiblePathogenIntegraionBkptPipelineSpark.class);

    @Argument(doc = "sam file for aligned contigs", shortName = "contigSAMFile",
            fullName = "contigSAMFile")
    private String outputAssemblyAlignments;

    @Argument(doc = "filename for output vcf", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String vcfOutputFileName;


    @ArgumentCollection
    private StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection evidenceAndAssemblyArgs
            = new StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection();

    @ArgumentCollection
    private StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection discoverStageArgs
            = new StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection();

    @Argument(doc = "size of the maximum gap (uncovered by alignments) a read must have for it to be included in output", shortName = "mgs",
            fullName = "maxGapSize")
    private int maxGapSize;

    @Argument(doc = "whether to write alignments of all successful assemblies", shortName = "allAln",
            fullName = "allAln")
    private boolean writeAllAlignments = false;

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

        final SAMFileHeader headerForReads = getHeaderForReads();

        final JavaRDD<GATKRead> contigAlignments = gather(ctx);
        if (contigAlignments.isEmpty()) return;

        final JavaRDD<AlignedContig> selectedAlignments =
                ExtractSuspectedPathogenAlignmentsSpark.selectMappedContigs(contigAlignments, maxGapSize, headerForReads, localLogger);
        if(selectedAlignments.isEmpty()) {
            localLogger.warn("Could not find any sites suggesting pathogen integration site.");
            return;
        }

        // write sites to vcf
        final SAMSequenceDictionary referenceSequenceDictionary = headerForReads.getSequenceDictionary();
        final Broadcast<ReferenceMultiSource> referenceBroadcast = ctx.broadcast(getReference());
        inferIntegrationSitesAndWriteVCF(selectedAlignments, referenceBroadcast, referenceSequenceDictionary);
    }

    private JavaRDD<GATKRead> gather(final JavaSparkContext ctx) {
        final JavaRDD<GATKRead> reads = getUnfilteredReads();
        final SAMFileHeader headerForReads = getHeaderForReads();
        final String outputName = writeAllAlignments ? outputAssemblyAlignments.replace(".sam", ".all.sam").replace(".bam", ".all.bam") : null;

        final List<String> refNames = AlignedAssemblyOrExcuse.getRefNames(headerForReads);

        // gather evidence, run assembly, and align
        final List<GATKRead> contigAlignments =
                FindBreakpointEvidenceSpark.gatherEvidenceAndWriteContigSamFile(ctx, evidenceAndAssemblyArgs, headerForReads, reads, outputName, localLogger)
                        .getAlignedAssemblyOrExcuseList().stream()
                        .filter(AlignedAssemblyOrExcuse::isNotFailure)
                        .flatMap(aaoe -> aaoe.toSAMStreamForAlignmentsOfThisAssembly(headerForReads, refNames))
                        .map(SAMRecordToGATKReadAdapter::new)
                        .collect(Collectors.toList());
        return ctx.parallelize(contigAlignments);
    }

    private void inferIntegrationSitesAndWriteVCF(final JavaRDD<AlignedContig> selectedAlignments,
                                                  final Broadcast<ReferenceMultiSource> broadcastReference,
                                                  final SAMSequenceDictionary referenceSequenceDictionary) {
        final JavaRDD<VariantContext> annotatedVariants =
                selectedAlignments
                        .flatMapToPair(DetectPossiblePathogenIntegraionBkptPipelineSpark::turnToNARL)
                        .map(noveltyTypeAndEvidence ->
                                DiscoverVariantsFromContigAlignmentsSAMSpark.annotateVariant(
                                        noveltyTypeAndEvidence._1,
                                        null,
                                        noveltyTypeAndEvidence._2._1,
                                        noveltyTypeAndEvidence._2._2,
                                        broadcastReference));
        // as a intermediate tool, add annotation in output VCF for the number of alignments that evidence contigs have and signature (needs enum)

        SVVCFWriter.writeVCF(annotatedVariants.collect(), vcfOutputFileName, referenceSequenceDictionary, localLogger);
    }

    // TODO: 10/11/17 add a function to extract signature of source contigs alignment signature
    @SuppressWarnings("unchecked")
    private static Iterator<Tuple2<NovelAdjacencyReferenceLocations, Tuple2<SvType, Iterable<ChimericAlignment>>>> turnToNARL(final AlignedContig contig) {

        // for a contig, depending on if it has:
        // one: head and tail clipping leads to BND records with "inserted sequence" in ALT field
        // two: BND records {chr BND, SS BND, OrderSwitch BND--look for inserted sequence,
        // two: DEL records--must be large substitution
        // 2+ & and ambiguous configuration: only output alignments

        if (contig.hasEquallyGoodAlnConfigurations || SvDiscoverFromLocalAssemblyContigAlignmentsSpark.hasOnly2Alignments(contig)) {
            return EmptyListIterator.INSTANCE;
        }
        final List<Tuple2<NovelAdjacencyReferenceLocations, Tuple2<SvType, Iterable<ChimericAlignment>>>> result;
        if (contig.alignmentIntervals.size() == 1) {

        } else {
            if ( SvDiscoverFromLocalAssemblyContigAlignmentsSpark.isSameChromosomeMapping(contig) ) { // same chr
                if (SvDiscoverFromLocalAssemblyContigAlignmentsSpark.isLikelyInvBreakpointOrInsInv(contig)) { // SS

                } else if (SvDiscoverFromLocalAssemblyContigAlignmentsSpark.isSuggestingRefBlockOrderSwitch(contig)){ // OS

                } else { // ins del
                    return new SingletonListIterator(new Tuple2<>(contig.contigSequence, ChimericAlignment.parseOneContig(contig, DEFAULT_MIN_ALIGNMENT_LENGTH)));
                }
            } else { // diff chr

            }
        }

        return null;
    }
}
