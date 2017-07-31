package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.denoising.svd;

import org.apache.commons.math3.linear.RealMatrix;

/**
 * TODO
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class RandomizedSingularValueDecomposer {
    private static final int RANDOM_SEED = 1216;

    private final RealMatrix matrix;
    private final int numEigenvalues;
    private final int numOversamples;
    private final int numPowerIterations;
    private final PowerIterationNormalizer.Mode powerIterationNormalizerMode;
    private final boolean transpose;
    private final boolean flipSign;

    private static final int NUM_OVERSAMPLES_DEFAULT = 10;
    private static final double SMALL_NUM_EIGENVALUES_FACTOR = 0.1;
    private static final int NUM_POWER_ITERATIONS_DEFAULT = 4;
    private static final int NUM_POWER_ITERATIONS_DEFAULT_FOR_SMALL_NUM_EIGENVALUES = 7;
    private static final boolean FLIP_SIGN_DEFAULT = true;

    public RandomizedSingularValueDecomposer(final RealMatrix matrix,
                                             final int numEigenvalues) {
        this(matrix, numEigenvalues,
                NUM_OVERSAMPLES_DEFAULT,
                numEigenvalues < SMALL_NUM_EIGENVALUES_FACTOR * Math.min(matrix.getRowDimension(), matrix.getColumnDimension()) ?
                        NUM_POWER_ITERATIONS_DEFAULT_FOR_SMALL_NUM_EIGENVALUES :
                        NUM_POWER_ITERATIONS_DEFAULT,
                PowerIterationNormalizer.Mode.AUTO,
                matrix.getRowDimension() < matrix.getColumnDimension(),
                FLIP_SIGN_DEFAULT);
    }

    public RandomizedSingularValueDecomposer(final RealMatrix matrix,
                                             final int numEigenvalues,
                                             final int numOversamples,
                                             final int numPowerIterations,
                                             final PowerIterationNormalizer.Mode powerIterationNormalizerMode,
                                             final boolean transpose,
                                             final boolean flipSign) {
        this.matrix = matrix;
        this.numEigenvalues = numEigenvalues;
        this.numOversamples = numOversamples;
        this.numPowerIterations = numPowerIterations;
        this.powerIterationNormalizerMode = powerIterationNormalizerMode;
        this.transpose = transpose;
        this.flipSign = flipSign;
    }

    private static final class PowerIterationNormalizer {
        private enum Mode {
            AUTO, QR, LU, NONE
        }
    }

}
