package org.broadinstitute.hellbender.tools.copynumber.gcnv;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.Utils;

import javax.annotation.Nullable;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

/**
 * Helper class for single sample gCNV postprocessing
 */
public class GCNVPostProcessor {

    /**
     * VCF miscellaneous header keys
     */
    public static final String SAMPLE_NAME = "sample_name";

    /**
     * VCF FORMAT header keys
     */
    public static final String CN_MAP = "CNMAP";
    public static final String CN_MEAN = "CNMEAN";
    public static final String CN_STD = "CNSTD";
    public static final String CNQ = "CNQ";

    /**
     * Cap all quality scores to this value
     */
    public static final int MAX_QUAL_SCORE = Integer.MAX_VALUE;

    /**
     * Precision of Phred-scale scores
     */
    public static final double PHRED_SCORE_PRECISION = 0.1;

    /**
     * Maximum tolerated log probability value. Values between this constant and 0.0 as considered as a 0.0.
     * Values above this threshold are considered to indicate a log probability calculation problem that is
     * worth to be reported as an error or warning.
     */
    public static final double MAX_LOG_PROB = 1e-3;

    private final List<Allele> alleles;
    private final IntegerCopyNumberStateCollection integerCopyNumberStateCollection;
    private final String sampleName;
    private final VariantContextWriter outputWriter;

    protected GCNVPostProcessor(final VariantContextWriter outputWriter,
                      final IntegerCopyNumberStateCollection integerCopyNumberStateCollection,
                      final String sampleName) {
        this.outputWriter = Utils.nonNull(outputWriter);
        this.integerCopyNumberStateCollection = Utils.nonNull(integerCopyNumberStateCollection);
        this.alleles = integerCopyNumberStateCollection.getAlleles();
        this.sampleName = sampleName;
    }

    /**
     * Compose the header of the variant context
     *
     * @param commandLine command line of the VCF generation tool
     */
    protected void composeVariantContextHeader(@Nullable final String commandLine) {
        //TODO pass a list of intervals and add progress meter
        //ProgressMeter progressMeter = new ProgressMeter(1.0);
        //progressMeter.start();
        final VCFHeader result = new VCFHeader(Collections.emptySet(), Arrays.asList(sampleName));

        /* add VCF version */
        result.addMetaDataLine(new VCFHeaderLine(VCFHeaderVersion.VCF4_2.getFormatString(),
                VCFHeaderVersion.VCF4_2.getVersionString()));

        /* add command line */
        if (commandLine != null) {
            result.addMetaDataLine(new VCFHeaderLine("command", commandLine));
        }

        /* add miscellaneous header keys */
        result.addMetaDataLine(new VCFHeaderLine(SAMPLE_NAME, sampleName));

        /* header lines related to genotype formatting */
        result.addMetaDataLine(new VCFFormatHeaderLine(CN_MAP, 1,
                VCFHeaderLineType.Integer, "Copy number maximum a posteriori"));
        result.addMetaDataLine(new VCFFormatHeaderLine(CN_MEAN, 1,
                VCFHeaderLineType.Float, "Copy number posterior mean"));
        result.addMetaDataLine(new VCFFormatHeaderLine(CN_STD, 1,
                VCFHeaderLineType.Float, "Copy number posterior standard deviation"));
        result.addMetaDataLine(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_QUALITY_KEY, 1,
                VCFHeaderLineType.Integer, "Genotype call quality as the difference between the best and second best PL"));

        //TODO replace the VCFHeaderLineCount with A since we have a value per allele
        result.addMetaDataLine(new VCFFormatHeaderLine(CNQ, integerCopyNumberStateCollection.size(),
                VCFHeaderLineType.Float, "Posterior copy number vector"));

        /* INFO header lines*/
        result.addMetaDataLine(new VCFInfoHeaderLine(VCFConstants.END_KEY, 1,
                VCFHeaderLineType.Integer, "End coordinate of this variant"));
        outputWriter.writeHeader(result);
    }

    protected void writeChunkedVariantContext(final List<CopyNumberPosteriorLocatableRecord> copyNumberPosteriorLocatableRecordsChunk,
                                           final String variantPrefix) {
        for (CopyNumberPosteriorLocatableRecord posteriorRecord: copyNumberPosteriorLocatableRecordsChunk) {
            final VariantContext variantContext = composeVariantContext(posteriorRecord, variantPrefix);
            outputWriter.add(variantContext);
        }
    }

    /**
     *
     * @param copyNumberPosteriorLocatableRecord a posterior record to genotype
     * @param variantPrefix a variant prefix
     * @return composed variant context
     */
    protected VariantContext composeVariantContext(final CopyNumberPosteriorLocatableRecord copyNumberPosteriorLocatableRecord,
                                                 final String variantPrefix) {
        final VariantContextBuilder variantContextBuilder = new VariantContextBuilder();
        variantContextBuilder.alleles(alleles);
        variantContextBuilder.chr(copyNumberPosteriorLocatableRecord.getContig());
        variantContextBuilder.start(copyNumberPosteriorLocatableRecord.getStart());
        variantContextBuilder.stop(copyNumberPosteriorLocatableRecord.getEnd());
        variantContextBuilder.id(String.format(variantPrefix + "_%s_%d_%d",
                copyNumberPosteriorLocatableRecord.getContig(),
                copyNumberPosteriorLocatableRecord.getStart(),
                copyNumberPosteriorLocatableRecord.getEnd()));
        final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(sampleName);
        final List<Double> posteriorVector = getCopyNumberPosteriorVector(copyNumberPosteriorLocatableRecord);
        final int copyNumberMAP = calculateMAPCopyNumber(copyNumberPosteriorLocatableRecord);
        final double copyNumberMean = calculateMAPMeanCopyNumber(copyNumberPosteriorLocatableRecord);
        final double copyNumberStd = calculateMAPStandardDeviation(copyNumberPosteriorLocatableRecord, copyNumberMean);
        final int GQ = calculateGenotypeQuality(copyNumberPosteriorLocatableRecord);

        genotypeBuilder.attribute(CNQ, posteriorVector);
        genotypeBuilder.attribute(CN_MAP, copyNumberMAP);
        genotypeBuilder.attribute(CN_MEAN, copyNumberMean);
        genotypeBuilder.attribute(CN_STD, copyNumberStd);
        genotypeBuilder.GQ(GQ);
        //TODO add posterior copy number vector values as an attribute
        final Genotype genotype = genotypeBuilder.make();

        //Add allele information to the variant context
        variantContextBuilder.attribute(VCFConstants.END_KEY, copyNumberPosteriorLocatableRecord.getEnd());

        variantContextBuilder.genotypes(genotype);
        return variantContextBuilder.make();
    }

    @VisibleForTesting
    public List<Double> getCopyNumberPosteriorVector(final CopyNumberPosteriorLocatableRecord copyNumberPosteriorLocatableRecord) {
        return integerCopyNumberStateCollection.getCopyNumberStates().stream().mapToDouble(state ->
                copyNumberPosteriorLocatableRecord.getCopyNumberPosteriorFromLogScale(state))
                .boxed().collect(Collectors.toList());
    }

    @VisibleForTesting
    public int calculateMAPCopyNumber(final CopyNumberPosteriorLocatableRecord copyNumberPosteriorLocatableRecord) {
        final Optional<IntegerCopyNumberState> copyNumberStateMAP = integerCopyNumberStateCollection.getCopyNumberStates()
                .stream().max((a1, a2) -> Double.compare(copyNumberPosteriorLocatableRecord.getCopyNumberPosteriors().getCopyNumberPosterior(a1),
                copyNumberPosteriorLocatableRecord.getCopyNumberPosteriors().getCopyNumberPosterior(a2)));
        return copyNumberStateMAP.get().getCopyNumber();
    }

    @VisibleForTesting
    public double calculateMAPMeanCopyNumber(final CopyNumberPosteriorLocatableRecord copyNumberPosteriorLocatableRecord) {
        return integerCopyNumberStateCollection.getCopyNumberStates()
                .stream().mapToDouble(state -> state.getCopyNumber() *
                        copyNumberPosteriorLocatableRecord.getCopyNumberPosteriorFromLogScale(state)).sum();
    }

    @VisibleForTesting
    public double calculateMAPStandardDeviation(final CopyNumberPosteriorLocatableRecord copyNumberPosteriorLocatableRecord,
                                                 final double meanMAP) {
        final double variance = integerCopyNumberStateCollection.getCopyNumberStates().stream()
                .mapToDouble(state -> copyNumberPosteriorLocatableRecord.getCopyNumberPosteriorFromLogScale(state) *
                        FastMath.pow(meanMAP - state.getCopyNumber(), 2)).sum();
        return FastMath.sqrt(variance);
    }

    @VisibleForTesting
    public int calculateGenotypeQuality(final CopyNumberPosteriorLocatableRecord copyNumberPosteriorLocatableRecord) {
        final double[] sortedPosteriors = integerCopyNumberStateCollection.getCopyNumberStates().stream()
                .mapToDouble(state -> copyNumberPosteriorLocatableRecord.getCopyNumberPosteriors().getCopyNumberPosterior(state))
                .sorted().toArray();
        if (sortedPosteriors.length < 2) {
            throw new UserException.BadInput("There must be at least two copy number states included in the posterior output file");
        }
        final double bestQuality = GATKProtectedMathUtils.logProbToPhredScore(sortedPosteriors[sortedPosteriors.length - 1],
                true, MAX_QUAL_SCORE, MAX_LOG_PROB, PHRED_SCORE_PRECISION);
        final double secondBestQuality = GATKProtectedMathUtils.logProbToPhredScore(sortedPosteriors[sortedPosteriors.length - 2],
                true, MAX_QUAL_SCORE, MAX_LOG_PROB, PHRED_SCORE_PRECISION);
        return (int) (bestQuality - secondBestQuality);
    }

}
