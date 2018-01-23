package org.broadinstitute.hellbender.tools.copynumber;

import org.broadinstitute.hellbender.CommandLineProgramTest;

import java.io.File;

import static org.testng.Assert.*;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class ModelSegmentsIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_SUB_DIR = new File(toolsTestDir, "copynumber");
    private static final File TUMOR_DENOISED_COPY_RATIOS_FILE = new File(TEST_SUB_DIR, "model-segments-wes-tumor-denoised-copy-ratios-SM-74P4M-v1-chr20-downsampled.deduplicated.denoisedCR.tsv");
    private static final File TUMOR_ALLELIC_COUNTS_FILE = new File(TEST_SUB_DIR, "model-segments-wes-tumor-allelic-counts-SM-74P4M-v1-chr20-downsampled.deduplicated.allelicCounts.tsv");
    private static final File NORMAL_ALLELIC_COUNTS_FILE = new File(TEST_SUB_DIR, "model-segments-wes-normal-allelic-counts-SM-74NEG-v1-chr20-downsampled.deduplicated.allelicCounts.tsv");
}