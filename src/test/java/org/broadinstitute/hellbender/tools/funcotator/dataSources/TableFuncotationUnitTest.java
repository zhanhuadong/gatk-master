package org.broadinstitute.hellbender.tools.funcotator.dataSources;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.stream.Collectors;

/**
 * A Unit Test class for {@link TableFuncotation}
 * Created by jonn on 11/28/17.
 */
public class TableFuncotationUnitTest extends GATKBaseTest {

    //==================================================================================================================
    // Private Static Members:

    //==================================================================================================================
    // Private Members:

    //==================================================================================================================
    // Helper Methods:

    //==================================================================================================================
    // Data Providers:

    @DataProvider
    Object[][] provideForTestGet() {
        return new Object[][]{
                { "X" },
                { "A" },
                { "B" },
                { "Y" },
                { "L" },
                { "R" },
                { "Select" },
                { "Start" },
        };
    }

    @DataProvider
    Object[][] provideDataForTestSetFieldSerializationOverrideValue() {

        final TableFuncotation funcotation =
                new TableFuncotation(
                    Arrays.asList("A", "B", "C"), Arrays.asList("1", "2", "3"), "TableFuncotation"
                );

        return new Object[][] {
                { funcotation, "A", null },
                { funcotation, "A", "A" },
                { funcotation, "A", "B" },
                { funcotation, "A", "C" },
                { funcotation, "B", "A" },
                { funcotation, "B", "B" },
                { funcotation, "B", "C" },
                { funcotation, "C", "A" },
                { funcotation, "C", "B" },
                { funcotation, "C", "C" },
        };
    }

    @DataProvider
    Object[][] provideForTestSerializeToVcfString() {
        return new Object[][] {
                {
                    new TableFuncotation(Collections.emptyList(), Collections.emptyList(), "Empty"),
                    ""
                },
                {
                    new TableFuncotation(Collections.singletonList("A"), Collections.singletonList(("1")), "OneVal"),
                    "1"
                },
                {
                    new TableFuncotation(Arrays.asList("A", "B"), Arrays.asList("1", "2"), "TwoVals"),
                    "1|2"
                },
                {
                    new TableFuncotation(Arrays.asList("A", "B", "C"), Arrays.asList("1", "2", "3"), "ThreeVals"),
                    "1|2|3"
                },
        };
    }

    @DataProvider
    Object[][] provideListOfStrings() {
        return new Object[][] {
                { Collections.singletonList("A") },
                { Arrays.asList("A", "B") },
                { Arrays.asList("A", "B", "C") },
                { Arrays.asList("A", "B", "C", "D") },
                { Arrays.asList("A", "B", "C", "D", "E") },
        };
    }

    @DataProvider
    Object[][] provideForTestGetFieldNames() {
        //final TableFuncotation tableFuncotation, final LinkedHashSet<String> expected
        return new Object[][] {
                {
                    new TableFuncotation(Collections.emptyList(), Collections.emptyList(), "Empty"),
                    new LinkedHashSet<>(Collections.emptyList())
                },
                {
                    new TableFuncotation(Collections.singletonList("TESTFIELD"), Collections.singletonList("TESTVAL"), "OneField"),
                    new LinkedHashSet<>(Collections.singletonList("TESTFIELD"))
                },
                {
                    new TableFuncotation(Arrays.asList("TESTFIELD1", "TESTFIELD2"), Arrays.asList("TESTVAL1", "TESTVAL2"), "TwoFields"),
                    new LinkedHashSet<>(Arrays.asList("TESTFIELD1", "TESTFIELD2"))
                },
        };
    }

    @DataProvider
    Object[][] provideForTestGetField() {
        //final TableFuncotation tableFuncotation, final String fieldName, final String expected
        return new Object[][] {
                {
                    new TableFuncotation(Collections.singletonList("TESTFIELD"), Collections.singletonList("TESTVAL"), "OneField"),
                    "TESTFIELD",
                    "TESTVAL"
                },
                {
                    new TableFuncotation(Arrays.asList("TESTFIELD1", "TESTFIELD2"), Arrays.asList("TESTVAL1", "TESTVAL2"), "TwoFields"),
                    "TESTFIELD1",
                    "TESTVAL1"
                },
                {
                    new TableFuncotation(Arrays.asList("TESTFIELD1", "TESTFIELD2"), Arrays.asList("TESTVAL1", "TESTVAL2"), "TwoFields"),
                    "TESTFIELD2",
                    "TESTVAL2"
                },
                {
                    new TableFuncotation(Arrays.asList("TESTFIELD1", "TESTFIELD2", "TESTFIELD3"), Arrays.asList("TESTVAL1", "TESTVAL2", "TESTVAL3"), "ThreeFields"),
                    "TESTFIELD1",
                    "TESTVAL1"
                },
                {
                    new TableFuncotation(Arrays.asList("TESTFIELD1", "TESTFIELD2", "TESTFIELD3"), Arrays.asList("TESTVAL1", "TESTVAL2", "TESTVAL3"), "ThreeFields"),
                    "TESTFIELD2",
                    "TESTVAL2"
                },
                {
                    new TableFuncotation(Arrays.asList("TESTFIELD1", "TESTFIELD2", "TESTFIELD3"), Arrays.asList("TESTVAL1", "TESTVAL2", "TESTVAL3"), "ThreeFields"),
                    "TESTFIELD3",
                    "TESTVAL3"
                },
        };
    }

    @DataProvider
    Object[][] provideForTestGetFieldFail() {
        //final TableFuncotation tableFuncotation, final String fieldName, final String expected
        return new Object[][] {
                {
                    new TableFuncotation(Collections.emptyList(), Collections.emptyList(), "Empty"),
                    "TESTFIELD_OMICRON"
                },
                {
                    new TableFuncotation(Collections.singletonList("TESTFIELD"), Collections.singletonList("TESTVAL"), "OneField"),
                    "testfield"
                },
                {
                    new TableFuncotation(Collections.singletonList("TESTFIELD"), Collections.singletonList("TESTVAL"), "OneField"),
                    "table_TESTFIELD"
                },
                {
                    new TableFuncotation(Arrays.asList("TESTFIELD1", "TESTFIELD2", "TESTFIELD3"), Arrays.asList("TESTVAL1", "TESTVAL2", "TESTVAL3"), "ThreeFields"),
                    "TESTFIELD4"
                },
        };
    }

    //==================================================================================================================
    // Tests:

    @Test(dataProvider = "provideForTestGet")
    public void testGet(final String fieldValue) {

        final String fieldName = "PLACEHOLDER";

        final TableFuncotation funcotation = new TableFuncotation( Collections.singletonList(fieldName), Collections.singletonList(fieldValue), fieldName );
        Assert.assertEquals( funcotation.get(fieldName), fieldValue );
    }

    @Test(dataProvider = "provideDataForTestSetFieldSerializationOverrideValue")
    public void testSetFieldSerializationOverrideValue(final TableFuncotation funcotation,
                                                       final String fieldName,
                                                       String overrideValue) {
        if ( overrideValue == null ) {
            overrideValue = funcotation.get(fieldName);
        }
        else {
            funcotation.setFieldSerializationOverrideValue(fieldName, overrideValue);
        }
        Assert.assertEquals( funcotation.get(fieldName), overrideValue );
    }

    @Test(dataProvider = "provideForTestSerializeToVcfString")
    public void testSerializeToVcfString(final TableFuncotation funcotation,
                                         final String expected) {
        Assert.assertEquals( funcotation.serializeToVcfString(), expected );
    }

    @Test(dataProvider = "provideListOfStrings")
    public void testKeySet(final List<String> kvNames) {
        final TableFuncotation funcotation = new TableFuncotation(kvNames, kvNames.stream().map(s -> s + "VVV").collect(Collectors.toList()), "ListFuncotation");
        Assert.assertEquals( funcotation.keySet(), kvNames);
    }

    @Test(dataProvider = "provideListOfStrings")
    public void testValues(final List<String> kvNames) {
        final TableFuncotation funcotation = new TableFuncotation(kvNames.stream().map(s -> s + "KEY").collect(Collectors.toList()), kvNames, "ListFuncotation");
        Assert.assertEquals( funcotation.values(), kvNames);
    }

    @Test(dataProvider = "provideForTestGetFieldNames")
    public void testGetFieldNames(final TableFuncotation tableFuncotation, final LinkedHashSet<String> expected) {
        Assert.assertEquals(tableFuncotation.getFieldNames(), expected);
    }

    @Test(dataProvider = "provideForTestGetField")
    public void testGetField(final TableFuncotation tableFuncotation, final String fieldName, final String expected) {
        Assert.assertEquals(tableFuncotation.getField(fieldName), expected);
    }

    @Test(dataProvider = "provideForTestGetFieldFail", expectedExceptions = GATKException.class)
    public void testGetFieldFail(final TableFuncotation tableFuncotation, final String fieldName) {
        tableFuncotation.getField(fieldName);
    }

}
