package org.broadinstitute.hellbender.tools.copynumber.gcnv;

import htsjdk.samtools.util.Locatable;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Map;

/**
 * A single copy number posterior record for a specific interval
 */
public class CopyNumberPosteriorLocatableRecord implements Locatable {
    private final SimpleInterval interval;
    private final CopyNumberPosteriorRecord copyNumberPosteriorRecord;

    public CopyNumberPosteriorLocatableRecord(final SimpleInterval interval, final CopyNumberPosteriorRecord copyNumberStatePosteriors) {
        this.interval = Utils.nonNull(interval);
        this.copyNumberPosteriorRecord = Utils.nonNull(copyNumberStatePosteriors);
    }

    /**
     * Get the copy number posteriors for this record
     */
    public CopyNumberPosteriorRecord getCopyNumberPosteriors() {
        return copyNumberPosteriorRecord;
    }

    /**
     * TODO cache these results
     * Convert posterior probability back from log scale
     *
     * @param copyNumberState copy number state to calculate posterior probability for
     * @return posterior probability
     */
    public double getCopyNumberPosteriorFromLogScale(final IntegerCopyNumberState copyNumberState) {
        return FastMath.exp(this.getCopyNumberPosteriors().getCopyNumberPosterior(copyNumberState));
    }

    @Override
    public String getContig() {
        return interval.getContig();
    }

    @Override
    public int getStart() {
        return interval.getStart();
    }

    @Override
    public int getEnd() {
        return interval.getEnd();
    }
}
