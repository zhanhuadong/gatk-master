package org.broadinstitute.hellbender.tools.copynumber.gcnv;

import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import static org.testng.Assert.*;

/**
 * Unit tests for {@link GCNVPostProcessor}
 */
public class GCNVPostProcessorUnitTest {

    @DataProvider(name = "examplePosteriorRecords")
    public Object[][] testData() {
        final SimpleInterval interval = new SimpleInterval("1", 1, 10000);
        final List<IntegerCopyNumberState> copyNumberStateList = new ArrayList<>();
        copyNumberStateList.add(new IntegerCopyNumberState(0));
        copyNumberStateList.add(new IntegerCopyNumberState(1));
        copyNumberStateList.add(new IntegerCopyNumberState(2));
        copyNumberStateList.add(new IntegerCopyNumberState(3));
        final IntegerCopyNumberStateCollection copyNumberStateCollection =
                new IntegerCopyNumberStateCollection(copyNumberStateList.stream().map(s -> s.toString()).collect(Collectors.toList()));
        final CopyNumberPosteriorRecord posteriorRecord = new CopyNumberPosteriorRecord();
        final CopyNumberPosteriorLocatableRecord locatableRecord = new CopyNumberPosteriorLocatableRecord();


        return new Object[][] {
            {posteriorRecord}
        };
    }


    @Test(dataProvider = "examplePosteriorRecords")
    public void test(final CopyNumberPosteriorLocatableRecord posteriorLocatableRecord) {

    }
}