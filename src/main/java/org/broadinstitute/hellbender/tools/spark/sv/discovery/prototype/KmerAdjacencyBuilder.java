package org.broadinstitute.hellbender.tools.spark.sv.discovery.prototype;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.tools.spark.sv.utils.*;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVKmer.Base;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVFastqUtils.FastqRead;

import java.util.*;

@BetaFeature
@CommandLineProgramProperties(
        oneLineSummary = "(Internal) junk", summary = "complete crap",
        programGroup = StructuralVariantDiscoveryProgramGroup.class)
public class KmerAdjacencyBuilder extends CommandLineProgram {
    private static final long serialVersionUID = 1L;

    @Argument(doc = "input fastq",
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.INPUT_LONG_NAME)
    public String fastqFile;

    @Override protected Object doWork() {
        final List<FastqRead> reads = SVFastqUtils.readFastqFile(fastqFile);
        final int kSize = 63;
        final int minQ = 3;
        final int minKCount = 4;
        final int nKmers = reads.stream().mapToInt(read -> Math.min(0, read.getBases().length - kSize + 1)).sum();
        final Map<SVKmerLong, int[]> kmerCounts = new HashMap<>(SVUtils.hashMapCapacity(nKmers));
        for ( final FastqRead read : reads ) {
            final byte[] quals = read.getQuals();
            int trimLen;
            for ( trimLen = 0; trimLen != quals.length; ++trimLen ) {
                if ( quals[trimLen] < minQ ) break;
            }
            SVKmerizer.canonicalStream(Arrays.copyOf(read.getBases(),trimLen),kSize,new SVKmerLong(kSize))
                    .forEach(kmer -> {
                        final int[] counts = kmerCounts.computeIfAbsent((SVKmerLong)kmer, key -> new int[1]);
                        counts[0] += 1;
                    });
        }
        kmerCounts.entrySet().removeIf(entry -> entry.getValue()[0] < minKCount);
        final Map<SVKmerLong, int[]> kmerAdjacencyMap = new HashMap<>(SVUtils.hashMapCapacity(kmerCounts.size()));
        kmerCounts.forEach((kmer, value) -> {
            final SVKmerLong kkk = kmer.removeFirstAndLastBase(kSize);
            final int[] counts = kmerAdjacencyMap.computeIfAbsent(kkk, key -> new int[8]);
            counts[kmer.firstBase(kSize).ordinal()] += value[0];
            counts[kmer.lastBase().ordinal() + 4] += value[0];
        });
        for ( final Map.Entry<SVKmerLong, int[]> entry : kmerAdjacencyMap.entrySet() ) {
            final StringBuilder sb = new StringBuilder(entry.getKey().toString(kSize-2));
            Arrays.stream(entry.getValue()).forEach(iii -> {sb.append('\t'); sb.append(iii);});
            System.out.println(sb);
        }
        kmerAdjacencyMap.forEach((kmer, counts) -> {
            for (final Base base : Base.values()) {
                if (counts[base.ordinal()] > 0) {
                    final SVKmerLong predecessor = kmer.predecessor(base, kSize - 2).canonical(kSize - 2);
                    if (!kmerAdjacencyMap.containsKey(predecessor)) {
                        System.out.println("ZSource: " + predecessor.toString(kSize - 2));
                    }
                }
                if (counts[base.ordinal() + 4] > 0) {
                    final SVKmerLong successor = kmer.successor(base, kSize - 2).canonical(kSize - 2);
                    if (!kmerAdjacencyMap.containsKey(successor)) {
                        System.out.println("ZSink:   " + successor.toString(kSize - 2));
                    }
                }
            }
        });
        return null;
    }
}
