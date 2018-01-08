package org.broadinstitute.hellbender.engine.spark;

import org.broadinstitute.hellbender.engine.MultiIntervalShard;
import org.broadinstitute.hellbender.engine.Shard;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.Collections;
import java.util.Iterator;
import java.util.List;

public class ShardToMultiIntervalShardAdapter<T> implements MultiIntervalShard<T> {
    private final Shard<T> shard;

    public ShardToMultiIntervalShardAdapter(final Shard<T> shard) {
        this.shard = shard;
    }

    @Override
    public List<SimpleInterval> getIntervals() {
        return Collections.singletonList(shard.getInterval());
    }

    @Override
    public List<SimpleInterval> getPaddedIntervals() {
        return Collections.singletonList(shard.getPaddedInterval());
    }

    @Override
    public Iterator<T> iterator() {
        return shard.iterator();
    }
}
