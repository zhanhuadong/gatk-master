package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.engine.AssemblyRegionEvaluator;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.activityprofile.ActivityProfileState;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.realignmentfilter.Realigner;

import java.util.ArrayList;
import java.util.List;
import java.util.Optional;

public class SomaticActiveRegionEngine implements AssemblyRegionEvaluator {
    public static final int INDEL_START_QUAL = 30;
    public static final int INDEL_CONTINUATION_QUAL = 10;
    public static final double MAX_ALT_FRACTION_IN_NORMAL = 0.3;
    public static final int MAX_NORMAL_QUAL_SUM = 100;

    // if qual sum exceeds this amount, no need to continue realigning alt reads
    public static final int REALIGNMENT_QUAL_SUM_THRESHOLD = 80;

    public static final int MAX_REALIGNMENT_FAILS = 2;
    public static final int MINIMUM_BASE_QUALITY = 6;   // for active region determination
    public static final int MAX_MAPPING_QUALITY_TO_REALIGN = 55;

    final M2ArgumentCollection MTAC;
    final SAMFileHeader header;
    private final Optional<Realigner> realigner;
    private final boolean hasNormal;

    public SomaticActiveRegionEngine(final M2ArgumentCollection MTAC, final SAMFileHeader header) {
        this.MTAC = MTAC;
        hasNormal = MTAC.normalSampleName != null;
        this.header = header;
        realigner = MTAC.realignmentFilterArgumentCollection.bwaMemIndexImage == null ? Optional.empty() :
                Optional.of(new Realigner(MTAC.realignmentFilterArgumentCollection, header));
    }

    @Override
    public ActivityProfileState isActive(final AlignmentContext context, final ReferenceContext ref, final FeatureContext featureContext) {
        final byte refBase = ref.getBase();
        final SimpleInterval refInterval = ref.getInterval();

        if( context == null || context.getBasePileup().isEmpty() ) {
            return new ActivityProfileState(refInterval, 0.0);
        }

        final ReadPileup pileup = context.getBasePileup();
        final ReadPileup tumorPileup = pileup.getPileupForSample(MTAC.tumorSampleName, header);
        final List<Evidence> tumorAltReads = findPossibleAltReads(tumorPileup, refBase);
        final int tumorAltCount = tumorAltReads.size();
        final int tumorRefCount = tumorPileup.size() - tumorAltCount;
        final double tumorQualSum = tumorAltReads.stream().mapToDouble(ev -> ev.qual).sum();

        final double tumorLog10Odds = -QualityUtils.qualToErrorProbLog10(tumorQualSum) +
                MathUtils.log10Factorial(tumorAltCount) + MathUtils.log10Factorial(tumorRefCount) - MathUtils.log10Factorial(tumorPileup.size() + 1);

        if (tumorLog10Odds < MTAC.initialTumorLodThreshold) {
            return new ActivityProfileState(refInterval, 0.0);
        } else if (hasNormal) {
            final ReadPileup normalPileup = pileup.getPileupForSample(MTAC.normalSampleName, header);
            final List<Evidence> normalAltReads = findPossibleAltReads(normalPileup, refBase);
            final int normalAltCount = normalAltReads.size();
            final double normalQualSum = normalAltReads.stream().mapToDouble(ev -> ev.qual).sum();
            if (normalAltCount > normalPileup.size() * MAX_ALT_FRACTION_IN_NORMAL && normalQualSum > MAX_NORMAL_QUAL_SUM) {
                return new ActivityProfileState(refInterval, 0.0);
            }
        } else {
            final List<VariantContext> germline = featureContext.getValues(MTAC.germlineResource, refInterval);
            if (!germline.isEmpty() && germline.get(0).getAttributeAsDoubleList(VCFConstants.ALLELE_FREQUENCY_KEY, 0.0).get(0) > MTAC.maxPopulationAlleleFrequency) {
                return new ActivityProfileState(refInterval, 0.0);
            }
        }

        if (!MTAC.genotypePonSites && !featureContext.getValues(MTAC.pon, new SimpleInterval(context.getContig(), (int) context.getPosition(), (int) context.getPosition())).isEmpty()) {
            return new ActivityProfileState(refInterval, 0.0);
        }

        if(realigner.isPresent() && calculateEvidenceAfterRealignment(pileup, tumorAltReads, tumorRefCount) < MTAC.initialTumorLodThreshold) {
            return new ActivityProfileState(refInterval, 0.0);
        }

        return new ActivityProfileState( refInterval, 1.0, ActivityProfileState.Type.NONE, null);
    }

    private double calculateEvidenceAfterRealignment(ReadPileup pileup, List<Evidence> tumorAltReads, int tumorRefCount) {
        final Locatable loc = pileup.getLocation();
        final Interval interval = new Interval(loc.getContig(), loc.getStart(), loc.getEnd());
        final Interval supposedRealignmentLocation = realigner.get().getRealignmentCoordinates(interval);

        int realignmentFailCount = 0;
        double qualSum = 0;
        int altCount = 0;

        for (final Evidence ev : tumorAltReads) {
            if (Realigner.isMarkedAsFailedRealignment(ev.read)) {
                realignmentFailCount++;
                continue;
            }

            if (ev.read.getMappingQuality() >= MAX_MAPPING_QUALITY_TO_REALIGN || qualSum > REALIGNMENT_QUAL_SUM_THRESHOLD || realigner.get().mapsToSupposedLocation(ev.read, supposedRealignmentLocation)) {
                qualSum += ev.qual;
                altCount++;
            } else {
                Realigner.markFailedRealignment(ev.read);
                realignmentFailCount++;
                if (realignmentFailCount > MAX_REALIGNMENT_FAILS) {
                    qualSum = 0;
                    break;
                }
            }
        }

        return -QualityUtils.qualToErrorProbLog10(qualSum) +
                MathUtils.log10Factorial(altCount) + MathUtils.log10Factorial(tumorRefCount) - MathUtils.log10Factorial(altCount + tumorRefCount + 1);
    }

    private static int getCurrentOrFollowingIndelLength(final PileupElement pe) {
        return pe.isDeletion() ? pe.getCurrentCigarElement().getLength() : pe.getLengthOfImmediatelyFollowingIndel();
    }

    private static double indelQual(final int indelLength) {
        return INDEL_START_QUAL + (indelLength - 1) * INDEL_CONTINUATION_QUAL;
    }

    private List<Evidence> findPossibleAltReads(final ReadPileup pileup, final byte refBase) {
        //TODO: replace with buffer?
        final List<Evidence> result = new ArrayList<>();
        for (final PileupElement pe : pileup) {
            final double altQual = altQuality(pe, refBase);
            if (altQual > 0) {
                result.add(new Evidence(pe.getRead(), altQual));
            }
        }

        return result;
    }

    private static double altQuality(final PileupElement pe, final byte refBase) {
        final int indelLength = getCurrentOrFollowingIndelLength(pe);
        if (indelLength > 0) {
            return indelQual(indelLength);
        } else if (isNextToUsefulSoftClip(pe)) {
            return indelQual(1);
        } else if (pe.getBase() != refBase && pe.getQual() > MINIMUM_BASE_QUALITY) {
            return pe.getQual();
        } else {
            return 0;
        }
    }

    // check that we're next to a soft clip that is not due to a read that got out of sync and ended in a bunch of BQ2's
    // we only need to check the next base's quality
    private static boolean isNextToUsefulSoftClip(final PileupElement pe) {
        final int offset = pe.getOffset();
        return pe.getQual() > MINIMUM_BASE_QUALITY &&
                ((pe.isBeforeSoftClip() && pe.getRead().getBaseQuality(offset + 1) > MINIMUM_BASE_QUALITY)
                        || (pe.isAfterSoftClip() && pe.getRead().getBaseQuality(offset - 1) > MINIMUM_BASE_QUALITY));
    }

    private static class Evidence {
        final public GATKRead read;
        final public double qual;

        public Evidence(final GATKRead read, final double qual) {
            this.read = read;
            this.qual = qual;
        }
    }

}
