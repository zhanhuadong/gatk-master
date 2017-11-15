package org.broadinstitute.hellbender.utils.realignmentfilter;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.liftover.LiftOver;
import htsjdk.samtools.util.Interval;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAligner;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

public class Realigner {
    public static final int MIN_MAP_QUALITY_FOR_REALIGNED_READS = 40;


    private final BwaMemAligner aligner;
    private final Optional<LiftOver> liftoverToRealignmentCoordinates;
    private final Function<Interval, Interval> convertToRealignmentCoordinates;
    private final List<String> realignmentContigs;

    public Realigner(final RealignmentFilterArgumentCollection rfac, final SAMFileHeader bamHeader) {
        final BwaMemIndex index = new BwaMemIndex(rfac.bwaMemIndexImage);
        realignmentContigs = index.getReferenceContigNames();
        aligner = new BwaMemAligner(index);
        liftoverToRealignmentCoordinates = rfac.liftoverChainFile == null ? Optional.empty() : Optional.of(new LiftOver(rfac.liftoverChainFile));
        convertToRealignmentCoordinates = liftoverToRealignmentCoordinates.isPresent() ? loc -> liftoverToRealignmentCoordinates.get().liftOver(loc) : loc -> loc;

        aligner.setMinSeedLengthOption(rfac.minSeedLength);
        aligner.setDropRatioOption((float) rfac.dropRatio);
        aligner.setSplitFactorOption((float) rfac.splitFactor);
    }

    //TODO: make a dedicated SAM read tag
    public static void markFailedRealignment(final GATKRead read) {
        read.setMappingQuality(0);
    }

    //TODO: same comment
    public static boolean isMarkedAsFailedRealignment(final GATKRead read) {
        return read.getMappingQuality() == 0;
    }

    public Interval getRealignemntCoordinates(final Interval bamCoordinates) {
        return convertToRealignmentCoordinates.apply(bamCoordinates);
    }


    public boolean mapsToSupposedLocation(final GATKRead read, final Interval supposedRealignmentLocation) {

        final List<BwaMemAlignment> alignments = aligner.alignSeqs(Arrays.asList(read), GATKRead::getBases).get(0);
        if (supposedRealignmentLocation == null || alignments.isEmpty()) {
            return false;
        }

        if (alignments.size() > 1) {
            // sort by descending mapping quality
            Collections.sort(alignments, (a1,a2) -> Double.compare(a2.getMapQual(), a1.getMapQual()));
        }

        final BwaMemAlignment alignment = alignments.get(0);
        if (alignment.getMapQual() < MIN_MAP_QUALITY_FOR_REALIGNED_READS) {
            return false;
        }

        final int contigId = alignment.getRefId();
        if (contigId < 0) {
            return false;
        }

        // TODO: is this necessary?   Maybe we could just check the contig
        // TODO: without checking all the coordinates.  We could also check the contig within some liftover fudge factor of 5 megabases or so
        // TODO: then we wouldn't need a liftover
        return new Interval(realignmentContigs.get(contigId), alignment.getRefStart(), alignment.getRefEnd()).overlaps(supposedRealignmentLocation);
    }
}
