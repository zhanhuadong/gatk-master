package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.List;

/**
 * An interface to represent shards of arbitrary data spanning multiple intervals.
 *
 * @param <T> Type of data in this shard
 */
public interface MultiIntervalShard<T> extends Iterable<T>, Shard<T>{

    /**
     * @return A List of this shard's intervals
     */
    List<SimpleInterval> getIntervals();

    /**
     * @return A List of this shard's intervals, with padding added to each interval on both sides
     */
    List<SimpleInterval> getPaddedIntervals();

    @Override
    default SimpleInterval getInterval(){
        return IntervalUtils.getSpanningInterval(getIntervals());
    }

    @Override
    default SimpleInterval getPaddedInterval(){
        return IntervalUtils.getSpanningInterval(getPaddedIntervals());
    }
}
