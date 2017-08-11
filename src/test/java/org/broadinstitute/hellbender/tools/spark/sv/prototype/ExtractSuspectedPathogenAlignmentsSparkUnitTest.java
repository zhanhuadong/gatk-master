package org.broadinstitute.hellbender.tools.spark.sv.prototype;

import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.prototype.AlnModType;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import scala.Tuple2;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class ExtractSuspectedPathogenAlignmentsSparkUnitTest extends BaseTest {

    @DataProvider(name = "testAlignmentsCoverage")
    private Object[][] CreateTestDataForestAlignmentsCoverage() {
        final List<Object[]> data = new ArrayList<>(20);

        AlignmentInterval one = new AlignmentInterval(new SimpleInterval("chr19", 108240, 108401), 1, 160, TextCigarCodec.decode("160M"), false, 0, 0, 160, AlnModType.NONE);
        AlignedContig tig = new AlignedContig("asm000000:tig00000", "GAATAATCTATGGCATGAAAGATTTTACCTGTCAACAGTGGCTGGCTCTTCATGGTTGCTACAATGAGTGTGTAAGATTCTGAAGAACTCCTTTAATAAGCCTAAACTTAATGTTCAACTTAGAATAAATACAATTCTTCTAAATTTTTTTGAATAATTT".getBytes(),
                Arrays.asList(one), false);
        data.add(new Object[]{tig, 160, 0, false});

        one = new AlignmentInterval(new SimpleInterval("chr5", 66018509, 66018581), 180, 248, TextCigarCodec.decode("179H35M4D34M7H"), true, 0, 5, 44, AlnModType.NONE);
        AlignmentInterval two = new AlignmentInterval(new SimpleInterval("chr1", 66379, 66452), 138, 211, TextCigarCodec.decode("137S74M44S"), false, 60, 0, 74, AlnModType.NONE);
        AlignmentInterval three = new AlignmentInterval(new SimpleInterval("chr1", 66373, 66483), 1, 116, TextCigarCodec.decode("41M5I70M139H"), false, 9, 11, 60, AlnModType.NONE);
        tig = new AlignedContig("asm000000:tig00012", "TATATATTATTATATAATATAATATATATTATATAATATATTTTATTATATAATATAATATATATTATATAATATAATATATTTTATTATATAAATATATATTATATTATATAATATAATATATATTTATATAATATAAATATATATTATATTATATAATATAATATATATTATATAATATATTTTATTATATAAATATATATTATATTATATAATATAATATATATTTTATTATATAATATATATTATATATTT".getBytes(),
                Arrays.asList(three, two, one), false);
        data.add(new Object[]{tig, 227, 21, false});

        one = new AlignmentInterval(new SimpleInterval("chrUn_KN707875v1_decoy", 1835, 1945), 74, 184, TextCigarCodec.decode("73H111M376H"), true, 3, 11, 56, AlnModType.NONE);
        two = new AlignmentInterval(new SimpleInterval("chr14_GL000225v1_random", 64907, 65010), 81, 184, TextCigarCodec.decode("376H104M80H"), false, 22, 11, 49, AlnModType.NONE);
        three = new AlignmentInterval(new SimpleInterval("chrUn_JTFH01000971v1_decoy", 1493, 1568), 208, 283, TextCigarCodec.decode("207H76M277H"), true, 0, 8, 36, AlnModType.NONE);
        AlignmentInterval four = new AlignmentInterval(new SimpleInterval("chr18", 110284, 110484), 361, 560, TextCigarCodec.decode("141M1D59M360S"), false, 29, 16, 108, AlnModType.NONE);

        tig = new AlignedContig("asm031558:tig00078", "CATGTGCATTGATTCCATTCCATTCCTTTCAATTCCAAGTCTCCTGGCTGCCACTGCCACAGTGCAGAGGTCTGCTGAGATGCACGGGAGCCCGCCGTCCTCTCTCTGCCCATGTCCGTAAGTGAAGTTCTGGCCAGGGCTCCCCACGGTGGCTCTCCCGACACCTTCGGACGGCTCCCTCCCCCAGTAGTCCGCGGAAGGGAACAACGAAGGAGACTCGTTTGAACTCCGAGCCAAAGCGAAACCCTGCTAGCCTGCTTTGAGCAGAACATGTACCGGGGCAACACCACCAGAAAACAGCTGGCCAAGACCATTGGCATTTCGGATGCTAGGGTCCACATTTGGTTTCAGAATGAGAGGATCACCTAGGTGATCAGTGCAGTGATATGTCACACAAATTCCTATAGAGAGAGCCTAGAAAGCCTTACATCACCTGGGTGATCAGTGCAGATATCTGACACAATGCCCCAATATACAGTGCCTAGACAAGACTTCCATCACCTGGGTGATCAGTGAGGAGGTATATCACAAATCCCCCTCTAGGCAGAGAATAGAGAAGA".getBytes(),
                Arrays.asList(one, two, three, four), false);
        data.add(new Object[]{tig, 387, 176, true});

        return data.toArray(new Object[data.size()][]);
    }

    @Test(dataProvider = "testAlignmentsCoverage", groups = "sv")
    public void testAlignmentsCoverage(final AlignedContig input, final int expectedCoverage, final int expectedMaxGapSz, final boolean expectedIsSuspect) {
        final Tuple2<Integer, Integer> coverageAndMaxGapSz = ExtractSuspectedPathogenAlignmentsSpark.alignmentsCoverage(input);
        final int coverage = coverageAndMaxGapSz._1,
                  maxGapSz = coverageAndMaxGapSz._2;
        Assert.assertEquals(coverage, expectedCoverage);
        Assert.assertEquals(maxGapSz, expectedMaxGapSz);
        Assert.assertEquals(ExtractSuspectedPathogenAlignmentsSpark.keepContigForPathSeqUse(input, 150), expectedIsSuspect);
    }

}
