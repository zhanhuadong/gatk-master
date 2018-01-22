package org.broadinstitute.hellbender.tools.funcotator.mafOutput;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.funcotator.DataSourceFuncotationFactory;
import org.broadinstitute.hellbender.tools.funcotator.Funcotation;
import org.broadinstitute.hellbender.tools.funcotator.Funcotator;
import org.broadinstitute.hellbender.tools.funcotator.OutputRenderer;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.stream.Collectors;

/**
 * A Funcotator output renderer for writing to MAF files.
 * Complies with version 2.4 of the MAF standard found here:
 * https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4
 * Created by jonn on 12/5/17.
 */
public class MafOutputRenderer extends OutputRenderer {

    //==================================================================================================================
    // Public Static Members:

    /**
     * Version of the MAF standard that this {@link MafOutputRenderer} writes.
     */
    public static String VERSION = "2.4";

    //==================================================================================================================
    // Private Static Members:

    /**
     * The string representing a comment in a MAF file.
     */
    protected static final String COMMENT_STRING = "#";

    /**
     * Value to insert into columns with unspecified annotations.
     */
    protected static final String UNKNOWN_VALUE_STRING = "__UNKNOWN__";

    /**
     * Value to insert into unused annotation columns.
     */
    protected static final String UNUSED_STRING = "NA";

    /**
     * Default set of columns to include in this {@link MafOutputRenderer}.
     * Order of the columns is preserved by the {@link LinkedHashMap}, while still being able to access each field via
     * the associated key.
     */
    private static final LinkedHashMap<String, String> defaultMap = new LinkedHashMap<>();

    /**
     * Map for: outputFieldName -> dataSourceFieldName1, dataSourceFieldName2 ...
     *
     * This map informs how to fill out outputFields using the input funcotations.
     *
     * Strings in the value list are in order of priority for use in the outputField.
     * That is, if both dataSourceFieldName1 and dataSourceFieldName2 are present as funcotation fields, then only
     * dataSourceFieldName1 will be used as the outputField.
     */
    private static final Map<String, List<String>> outputFieldNameMap = new HashMap<>();

    /**
     * Delimiter for fields in the output MAF file.
     */
    private static final String FIELD_DELIMITER = "\t";

    //==================================================================================================================
    // Private Members:

    /**
     * {@link List} of the {@link DataSourceFuncotationFactory} objects that are being used in this run of {@link Funcotator}.
     */
    private final List<DataSourceFuncotationFactory> dataSourceFactories;

    /**
     * {@link Path} to the resulting MAF output file for this {@link MafOutputRenderer}.
     */
    private final Path outputFilePath;

    /**
     * {@link java.io.PrintWriter} to which to write the output MAF file.
     */
    private PrintWriter printWriter;

    //==================================================================================================================
    // Constructors:

    /**
     * Create a {@link MafOutputRenderer}.
     * @param outputFilePath {@link Path} to output file (must not be null).
     * @param dataSources {@link List} of {@link DataSourceFuncotationFactory} to back our annotations (must not be null).
     * @param unaccountedForDefaultAnnotations {@link LinkedHashMap} of default annotations that must be added.
     * @param unaccountedForOverrideAnnotations {@link LinkedHashMap} of override annotations that must be added.
     */
    public MafOutputRenderer(final Path outputFilePath,
                             final List<DataSourceFuncotationFactory> dataSources,
                             final LinkedHashMap<String, String> unaccountedForDefaultAnnotations,
                             final LinkedHashMap<String, String> unaccountedForOverrideAnnotations) {

        this.outputFilePath = outputFilePath;

        dataSourceFactories = dataSources;

        // Merge the annotations into our manualAnnotations:
        manualAnnotations = new LinkedHashMap<>();
        if ( unaccountedForDefaultAnnotations != null ) {
            manualAnnotations.putAll(unaccountedForDefaultAnnotations);
        }
        if ( unaccountedForOverrideAnnotations != null ) {
            manualAnnotations.putAll(unaccountedForOverrideAnnotations);
        }

        // Cache the manual annotation string so we can pass it easily into any Funcotations:
        manualAnnotationSerializedString = (manualAnnotations.size() != 0 ? String.join( FIELD_DELIMITER, manualAnnotations.values() ) + FIELD_DELIMITER : "");

        // Fill in our default output map:
        initializeDefaultMapWithKeys();

        // Fill in our default map for outputField -> funcotation name:
        initializeOutputFieldNameMap();
    }

    //==================================================================================================================
    // Override Methods:

    @Override
    public void open() {
        try {
            printWriter = new PrintWriter(Files.newOutputStream(outputFilePath));
            writeHeader();
        }
        catch (final IOException ex) {
            throw new UserException("Error opening output file path: " + outputFilePath.toUri().toString(), ex);
        }
    }

    @Override
    public void close() {
        printWriter.flush();
        printWriter.close();
    }

    @Override
    public void write(final VariantContext variant, final List<Funcotation> funcotations) {

        //TODO: For this to work correctly you'll need to group funcotations by reference / alternate allele pairs.

        // Make sure we only output the variant here if it passed all filters:
        if ( variant.isFiltered() ) {
            // We can ignore this variant since it was filtered out.
            return;
        }

        // Create our output maps:
        final LinkedHashMap<String, Object> outputMap = new LinkedHashMap<>(defaultMap);
        final LinkedHashMap<String, Object> extraFieldOutputMap = new LinkedHashMap<>();

        // Separate the funcotations into Gencode and Other:
        final List<GencodeFuncotation> gencodeFuncotations = new ArrayList<>(1);
        final List<Funcotation> otherFuncotations = new ArrayList<>(funcotations.size());
        for ( final Funcotation funcotation : funcotations ) {
            if ( GencodeFuncotation.class.isAssignableFrom(funcotation.getClass()) ) {
                gencodeFuncotations.add((GencodeFuncotation) funcotation);

                // We know we won't use secondary variant class here, so we add it to the extraFieldOutputMap:
                setField(extraFieldOutputMap, "secondary_variant_classification", ((GencodeFuncotation) funcotation).getSecondaryVariantClassification() );
            }
            else {
                // Add all the fields from the other funcotations into the extra field output:
                for ( final String field : funcotation.getFieldNames() ) {
                    setField(extraFieldOutputMap, field, funcotation.getField(field) );
                }
            }
        }

        // Set values for unused fields:
        outputMap.put("Score", UNUSED_STRING);
        outputMap.put("BAM_File", UNUSED_STRING);

        // Go through Gencode funcotations as they're the basis for each line:
        for ( final GencodeFuncotation gencodeFuncotation : gencodeFuncotations ) {
            setField(outputMap, "Hugo_Symbol", gencodeFuncotation.getHugoSymbol());
            setField(outputMap, "NCBI_Build", gencodeFuncotation.getNcbiBuild());
            setField(outputMap, "Chromosome", gencodeFuncotation.getChromosome());
            setField(outputMap, "Start_Position", gencodeFuncotation.getStart());
            setField(outputMap, "End_Position", gencodeFuncotation.getEnd());
            setField(outputMap, "Transcript_Strand", gencodeFuncotation.getTranscriptStrand());
            setField(outputMap, "VariantClassification", gencodeFuncotation.getVariantClassification());
            setField(outputMap, "Variant_Type", gencodeFuncotation.getVariantType());
            setField(outputMap, "Reference_Allele", gencodeFuncotation.getRefAllele());
            setField(outputMap, "Tumor_Seq_Allele1", gencodeFuncotation.getTumorSeqAllele1());
            setField(outputMap, "Tumor_Seq_Allele2", gencodeFuncotation.getTumorSeqAllele2());
            setField(outputMap, "Genome_Change", gencodeFuncotation.getGenomeChange() );
            setField(outputMap, "Annotation_Transcript", gencodeFuncotation.getAnnotationTranscript() );
            setField(outputMap, "Transcript_Exon", gencodeFuncotation.getTranscriptExonNumber() );
            setField(outputMap, "Transcript_Position", gencodeFuncotation.getTranscriptPos() );
            setField(outputMap, "cDNA_Change", gencodeFuncotation.getcDnaChange() );
            setField(outputMap, "Codon_Change", gencodeFuncotation.getCodonChange() );
            setField(outputMap, "Protein_Change", gencodeFuncotation.getProteinChange() );
            setField(outputMap, "gc_content", gencodeFuncotation.getGcContent() );
            setField(outputMap, "ref_context", gencodeFuncotation.getReferenceContext() );
            setField(outputMap, "Other_Transcripts", gencodeFuncotation.getOtherTranscripts() );
        }

        // Go through all output fields and see if any of the names in the value list are in our extraFieldOutputMap.
        // For any that match, we remove them from our extraFieldOutputMap and add them to the outputMap with the
        // correct key.
        for ( final Map.Entry<String, List<String>> entry : outputFieldNameMap.entrySet() ) {
            for ( final String fieldName : entry.getValue() ) {
                if ( extraFieldOutputMap.containsKey(fieldName) ) {
                    outputMap.put( entry.getKey(), extraFieldOutputMap.remove(fieldName) );
                }
            }
        }

        // Merge our output maps together:
        outputMap.putAll( extraFieldOutputMap );

        // Write the output:
        writeLine(
                outputMap.entrySet().stream()
                        .map(e -> e.getValue())
                        .map(Object::toString)
                        .collect(Collectors.joining("\t"))
        );

    }


    //==================================================================================================================
    // Static Methods:

    //==================================================================================================================
    // Instance Methods:

    /**
     * Set the field in the given {@code outputMap} specified by the given {@code key} to the given {@code value}.
     * Will only set this field if the given {@code value} is not null and when evaluated as a {@link String} is not empty.
     * @param outputMap {@link Map} of {@link String} to {@link Object} to hold output annotations.  Must not be {@code null}.
     * @param key {@link String} key to add to the output map.  Must not be {@code null}.
     * @param value {@link Object} value for the output map.
     */
    private void setField(final Map<String, Object> outputMap, final String key, final Object value) {
        Utils.nonNull(outputMap);
        Utils.nonNull(key);

        if ( (value != null) && (!value.toString().isEmpty())) {
            outputMap.put(key, value);
        }
        else {
            outputMap.put(key, UNKNOWN_VALUE_STRING);
        }
    }

    /**
     * Write the given line to the {@link #printWriter}.
     * @param line The {@link String} to write as a line to the {@link #printWriter}.
     */
    private void writeLine(final String line) {
        printWriter.write(line + System.lineSeparator());
    }

    /**
     * Write the header to the output file.
     */
    protected void writeHeader() {
        // Write out version:
        writeLine(COMMENT_STRING + "version " + VERSION);
        writeLine(COMMENT_STRING + COMMENT_STRING);

        // Write tool name and the data sources with versions:
        printWriter.write(COMMENT_STRING);
        printWriter.write(COMMENT_STRING);
        printWriter.write(" Funcotator ");
        printWriter.write(new SimpleDateFormat("yyyymmdd'T'hhmmss").format(new Date()));
        for (final DataSourceFuncotationFactory funcotationFactory : dataSourceFactories) {
            printWriter.write(" | ");
            printWriter.write(funcotationFactory.getName());
            printWriter.write(" ");
            printWriter.write(funcotationFactory.getVersion());
        }
        writeLine("");

        // Write the column headers:
        // TODO: Make sure we only account for the actual columns.  That is, we have a mapping of data sources to known columns then the other cols are appended at the end.  Make sure we don't double-count any of these!
        writeLine( defaultMap.keySet().stream().collect(Collectors.joining("\t")) );
    }

    /**
     * Initializes {@link MafOutputRenderer#defaultMap} with the default keys for the columns in a MAF file.
     */
    protected void initializeDefaultMapWithKeys() {
        
        // Baseline required fields:
        defaultMap.put("Hugo_Symbol",                                 UNKNOWN_VALUE_STRING );
        defaultMap.put("Entrez_Gene_Id",                              UNKNOWN_VALUE_STRING );
        defaultMap.put("Center",                                      UNKNOWN_VALUE_STRING );
        defaultMap.put("NCBI_Build",                                  UNKNOWN_VALUE_STRING );
        defaultMap.put("Chromosome",                                  UNKNOWN_VALUE_STRING );
        defaultMap.put("Start_Position",                              UNKNOWN_VALUE_STRING );
        defaultMap.put("End_Position",                                UNKNOWN_VALUE_STRING );
        defaultMap.put("Strand",                                      UNKNOWN_VALUE_STRING );
        defaultMap.put("Variant_Classification",                      UNKNOWN_VALUE_STRING );
        defaultMap.put("Variant_Type",                                UNKNOWN_VALUE_STRING );
        defaultMap.put("Reference_Allele",                            UNKNOWN_VALUE_STRING );
        defaultMap.put("Tumor_Seq_Allele1",                           UNKNOWN_VALUE_STRING );
        defaultMap.put("Tumor_Seq_Allele2",                           UNKNOWN_VALUE_STRING );
        defaultMap.put("dbSNP_RS",                                    UNKNOWN_VALUE_STRING );
        defaultMap.put("dbSNP_Val_Status",                            UNKNOWN_VALUE_STRING );
        defaultMap.put("Tumor_Sample_Barcode",                        UNKNOWN_VALUE_STRING );
        defaultMap.put("Matched_Norm_Sample_Barcode",                 UNKNOWN_VALUE_STRING );
        defaultMap.put("Match_Norm_Seq_Allele1",                      UNKNOWN_VALUE_STRING );
        defaultMap.put("Match_Norm_Seq_Allele2",                      UNKNOWN_VALUE_STRING );
        defaultMap.put("Tumor_Validation_Allele1",                    UNKNOWN_VALUE_STRING );
        defaultMap.put("Tumor_Validation_Allele2",                    UNKNOWN_VALUE_STRING );
        defaultMap.put("Match_Norm_Validation_Allele1",               UNKNOWN_VALUE_STRING );
        defaultMap.put("Match_Norm_Validation_Allele2",               UNKNOWN_VALUE_STRING );
        defaultMap.put("Verification_Status",                         UNKNOWN_VALUE_STRING );
        defaultMap.put("Validation_Status",                           UNKNOWN_VALUE_STRING );
        defaultMap.put("Mutation_Status",                             UNKNOWN_VALUE_STRING );
        defaultMap.put("Sequencing_Phase",                            UNKNOWN_VALUE_STRING );
        defaultMap.put("Sequence_Source",                             UNKNOWN_VALUE_STRING );
        defaultMap.put("Validation_Method",                           UNKNOWN_VALUE_STRING );
        defaultMap.put("Score",                                       UNKNOWN_VALUE_STRING );
        defaultMap.put("BAM_File",                                    UNKNOWN_VALUE_STRING );
        defaultMap.put("Sequencer",                                   UNKNOWN_VALUE_STRING );
        defaultMap.put("Tumor_Sample_UUID",                           UNKNOWN_VALUE_STRING );
        defaultMap.put("Matched_Norm_Sample_UUID",                    UNKNOWN_VALUE_STRING );

        // Required "optional" fields:
        defaultMap.put("Genome_Change",                               UNKNOWN_VALUE_STRING);
        defaultMap.put("Annotation_Transcript",                       UNKNOWN_VALUE_STRING);
        defaultMap.put("Transcript_Strand",                           UNKNOWN_VALUE_STRING);
        defaultMap.put("Transcript_Exon",                             UNKNOWN_VALUE_STRING);
        defaultMap.put("Transcript_Position",                         UNKNOWN_VALUE_STRING);
        defaultMap.put("cDNA_Change",                                 UNKNOWN_VALUE_STRING);
        defaultMap.put("Codon_Change",                                UNKNOWN_VALUE_STRING);
        defaultMap.put("Protein_Change",                              UNKNOWN_VALUE_STRING);
        defaultMap.put("Other_Transcripts",                           UNKNOWN_VALUE_STRING);
        defaultMap.put("Refseq_mRNA_Id",                              UNKNOWN_VALUE_STRING);
        defaultMap.put("Refseq_prot_Id",                              UNKNOWN_VALUE_STRING);
        defaultMap.put("SwissProt_acc_Id",                            UNKNOWN_VALUE_STRING);
        defaultMap.put("SwissProt_entry_Id",                          UNKNOWN_VALUE_STRING);
        defaultMap.put("Description",                                 UNKNOWN_VALUE_STRING);
        defaultMap.put("UniProt_AApos",                               UNKNOWN_VALUE_STRING);
        defaultMap.put("UniProt_Region",                              UNKNOWN_VALUE_STRING);
        defaultMap.put("UniProt_Site",                                UNKNOWN_VALUE_STRING);
        defaultMap.put("UniProt_Natural_Variations",                  UNKNOWN_VALUE_STRING);
        defaultMap.put("UniProt_Experimental_Info",                   UNKNOWN_VALUE_STRING);
        defaultMap.put("GO_Biological_Process",                       UNKNOWN_VALUE_STRING);
        defaultMap.put("GO_Cellular_Component",                       UNKNOWN_VALUE_STRING);
        defaultMap.put("GO_Molecular_Function",                       UNKNOWN_VALUE_STRING);
        defaultMap.put("COSMIC_overlap",                              UNKNOWN_VALUE_STRING);
        defaultMap.put("ping_mutations",                              UNKNOWN_VALUE_STRING);
        defaultMap.put("COSMIC_fusion_genes",                         UNKNOWN_VALUE_STRING);
        defaultMap.put("COSMIC_tissue_types_affected",                UNKNOWN_VALUE_STRING);
        defaultMap.put("COSMIC_total_alterations_in_gene",            UNKNOWN_VALUE_STRING);
        defaultMap.put("Tumorscape_Amplification_Peaks",              UNKNOWN_VALUE_STRING);
        defaultMap.put("Tumorscape_Deletion_Peaks",                   UNKNOWN_VALUE_STRING);
        defaultMap.put("TCGAscape_Amplification_Peaks",               UNKNOWN_VALUE_STRING);
        defaultMap.put("TCGAscape_Deletion_Peaks",                    UNKNOWN_VALUE_STRING);
        defaultMap.put("DrugBank",                                    UNKNOWN_VALUE_STRING);
        defaultMap.put("ref_context",                                 UNKNOWN_VALUE_STRING);
        defaultMap.put("gc_content",                                  UNKNOWN_VALUE_STRING);
        defaultMap.put("CCLE_ONCOMAP_overlapping_mutations",          UNKNOWN_VALUE_STRING);
        defaultMap.put("CCLE_ONCOMAP_total_mutations_in_gene",        UNKNOWN_VALUE_STRING);
        defaultMap.put("CGC_Mutation_Type",                           UNKNOWN_VALUE_STRING);
        defaultMap.put("CGC_Translocation_Partner",                   UNKNOWN_VALUE_STRING);
        defaultMap.put("CGC_Tumor_Types_Somatic",                     UNKNOWN_VALUE_STRING);
        defaultMap.put("CGC_Tumor_Types_Germline",                    UNKNOWN_VALUE_STRING);
        defaultMap.put("CGC_Other_Diseases",                          UNKNOWN_VALUE_STRING);
        defaultMap.put("DNARepairGenes_Role",                         UNKNOWN_VALUE_STRING);
        defaultMap.put("FamilialCancerDatabase_Syndromes",            UNKNOWN_VALUE_STRING);
        defaultMap.put("MUTSIG_Published_Results",                    UNKNOWN_VALUE_STRING);
        defaultMap.put("OREGANNO_ID",                                 UNKNOWN_VALUE_STRING);
        defaultMap.put("OREGANNO_Values",                             UNKNOWN_VALUE_STRING);
    }

    /**
     * Initializes the output map
     */
    private void initializeOutputFieldNameMap() {
        outputFieldNameMap.put( "Hugo_Symbol",                   Arrays.asList("Hugo_Symbol", "gene", "Gene") );
        outputFieldNameMap.put( "Entrez_Gene_Id",                Arrays.asList("Entrez_Gene_Id", "HGNC_Entrez Gene ID", "HGNC_Entrez Gene ID(supplied by NCBI)", "entrez_id", "gene_id") );
        outputFieldNameMap.put( "Center",                        Arrays.asList("Center", "center") );
        outputFieldNameMap.put( "NCBI_Build",                    Arrays.asList("NCBI_Build", "ncbi_build") );
        outputFieldNameMap.put( "Chromosome",                    Arrays.asList("Chromosome", "chr", "contig", "chromosome", "chrom", "Chrom") );
        outputFieldNameMap.put( "Start_position",                Arrays.asList("Start_position", "start", "Start", "start_pos", "pos") );
        outputFieldNameMap.put( "End_position",                  Arrays.asList("End_position", "end", "End", "end_pos") );
        outputFieldNameMap.put( "Strand",                        Arrays.asList("Strand", "Strand") );
        outputFieldNameMap.put( "Variant_Classification",        Arrays.asList("Variant_Classification", "variant_classification") );
        outputFieldNameMap.put( "Variant_Type",                  Arrays.asList("Variant_Type", "variant_type") );
        outputFieldNameMap.put( "Reference_Allele",              Arrays.asList("Reference_Allele", "ref", "ref_allele", "reference_allele") );
        outputFieldNameMap.put( "Tumor_Seq_Allele1",             Arrays.asList("Tumor_Seq_Allele1", "ref", "ref_allele", "reference_allele") );
        outputFieldNameMap.put( "Tumor_Seq_Allele2",             Arrays.asList("Tumor_Seq_Allele2", "alt", "alt_allele", "alt2", "alt_allele2", "alternate_allele2", "observed_allele2", "alternate_allele", "observed_allele", "alt1", "alt_allele1", "alternate_allele1", "observed_allele1") );
        outputFieldNameMap.put( "dbSNP_RS",                      Arrays.asList("dbSNP_RS", "dbsnp_rs") );
        outputFieldNameMap.put( "dbSNP_Val_Status",              Arrays.asList("dbSNP_Val_Status", "dbsnp_val_status") );
        outputFieldNameMap.put( "Tumor_Sample_Barcode",          Arrays.asList("Tumor_Sample_Barcode", "tumor_barcode", "tumor_id", "case_barcode", "case_id", "tumor_name") );
        outputFieldNameMap.put( "Matched_Norm_Sample_Barcode",   Arrays.asList("Matched_Norm_Sample_Barcode", "normal_barcode", "normal_id", "control_barcode", "control_id", "normal_name", "sample_name") );
        outputFieldNameMap.put( "Match_Norm_Seq_Allele1",        Arrays.asList("Match_Norm_Seq_Allele1", "Match_Norm_Seq_Allele1") );
        outputFieldNameMap.put( "Match_Norm_Seq_Allele2",        Arrays.asList("Match_Norm_Seq_Allele2", "Match_Norm_Seq_Allele2") );
        outputFieldNameMap.put( "Tumor_Validation_Allele1",      Arrays.asList("Tumor_Validation_Allele1", "Tumor_Validation_Allele1") );
        outputFieldNameMap.put( "Tumor_Validation_Allele2",      Arrays.asList("Tumor_Validation_Allele2", "Tumor_Validation_Allele2") );
        outputFieldNameMap.put( "Match_Norm_Validation_Allele1", Arrays.asList("Match_Norm_Validation_Allele1", "Match_Norm_Validation_Allele1") );
        outputFieldNameMap.put( "Match_Norm_Validation_Allele2", Arrays.asList("Match_Norm_Validation_Allele2", "Match_Norm_Validation_Allele2") );
        outputFieldNameMap.put( "Verification_Status",           Arrays.asList("Verification_Status", "Verification_Status") );
        outputFieldNameMap.put( "Validation_Status",             Arrays.asList("Validation_Status", "validation_status") );
        outputFieldNameMap.put( "Mutation_Status",               Arrays.asList("Mutation_Status", "status") );
        outputFieldNameMap.put( "Sequencing_Phase",              Arrays.asList("Sequencing_Phase", "phase") );
        outputFieldNameMap.put( "Sequence_Source",               Arrays.asList("Sequence_Source", "source") );
        outputFieldNameMap.put( "Validation_Method",             Arrays.asList("Validation_Method", "Validation_Method") );
        outputFieldNameMap.put( "Score",                         Arrays.asList("Score", "Score") );
        outputFieldNameMap.put( "BAM_file",                      Arrays.asList("BAM_file", "bam", "bam_file") );
        outputFieldNameMap.put( "Sequencer",                     Arrays.asList("Sequencer", "sequencer", "platform") );
        outputFieldNameMap.put( "Tumor_Sample_UUID",             Arrays.asList("Tumor_Sample_UUID", "tumor_uuid", "case_uuid", "tumor_barcode", "tumor_id", "case_barcode", "case_id", "tumor_name", "Tumor_Sample_Barcode") );
        outputFieldNameMap.put( "Matched_Norm_Sample_UUID",      Arrays.asList("Matched_Norm_Sample_UUID", "normal_uuid", "control_uuid", "normal_barcode", "normal_id", "control_barcode", "control_id", "normal_name", "sample_name", "Matched_Norm_Sample_Barcode") );
        outputFieldNameMap.put( "Genome_Change",                 Arrays.asList("Genome_Change", "genome_change") );
        outputFieldNameMap.put( "Annotation_Transcript",         Arrays.asList("Annotation_Transcript", "annotation_transcript", "transcript_id") );
        outputFieldNameMap.put( "Transcript_Strand",             Arrays.asList("Transcript_Strand", "transcript_strand") );
        outputFieldNameMap.put( "Transcript_Exon",               Arrays.asList("Transcript_Exon", "transcript_exon") );
        outputFieldNameMap.put( "Transcript_Position",           Arrays.asList("Transcript_Position", "transcript_position") );
        outputFieldNameMap.put( "cDNA_Change",                   Arrays.asList("cDNA_Change", "transcript_change") );
        outputFieldNameMap.put( "Codon_Change",                  Arrays.asList("Codon_Change", "codon_change") );
        outputFieldNameMap.put( "Protein_Change",                Arrays.asList("Protein_Change", "protein_change") );
        outputFieldNameMap.put( "Other_Transcripts",             Arrays.asList("Other_Transcripts", "other_transcripts") );
    }

    //==================================================================================================================
    // Helper Data Types:

    //------------------------------------------------------------------------------------------------------------------
    // Columns:
    
// Required:
// 1	Hugo_Symbol
// 2	Entrez_Gene_Id
// 3	Center
// 4	NCBI_Build
// 5	Chromosome
// 6	Start_Position
// 7	End_Position
// 8	Strand
// 9	Variant_Classification
//10	Variant_Type
//11	Reference_Allele
//12	Tumor_Seq_Allele1
//13	Tumor_Seq_Allele2
//14	dbSNP_RS
//15	dbSNP_Val_Status
//16	Tumor_Sample_Barcode
//17	Matched_Norm_Sample_Barcode
//18	Match_Norm_Seq_Allele1
//19	Match_Norm_Seq_Allele2
//20	Tumor_Validation_Allele1
//21	Tumor_Validation_Allele2
//22	Match_Norm_Validation_Allele1
//23	Match_Norm_Validation_Allele2
//24	Verification_Status
//25	Validation_Status
//26	Mutation_Status
//27	Sequencing_Phase
//28	Sequence_Source
//29	Validation_Method
//30	Score
//31	BAM_File
//32	Sequencer
//33    Tumor_Sample_UUID
//34    Matched_Norm_Sample_UUID

// "Optional" Required:
//35    Genome_Change
//36    Annotation_Transcript
//37    Transcript_Strand
//38    Transcript_Exon
//39    Transcript_Position
//40    cDNA_Change
//41    Codon_Change
//42    Protein_Change
//43    Other_Transcripts
//44    Refseq_mRNA_Id
//45    Refseq_prot_Id
//46    SwissProt_acc_Id
//47    SwissProt_entry_Id
//48    Description
//49    UniProt_AApos
//50    UniProt_Region
//51    UniProt_Site
//52    UniProt_Natural_Variations
//53    UniProt_Experimental_Info
//54    GO_Biological_Process
//55    GO_Cellular_Component
//56    GO_Molecular_Function
//57    COSMIC_overlap
//58    ping_mutations
//59    COSMIC_fusion_genes
//60    COSMIC_tissue_types_affected
//61    COSMIC_total_alterations_in_gene
//62    Tumorscape_Amplification_Peaks
//63    Tumorscape_Deletion_Peaks
//64    TCGAscape_Amplification_Peaks
//65    TCGAscape_Deletion_Peaks
//66    DrugBank
//67    ref_context
//68    gc_content
//69    CCLE_ONCOMAP_overlapping_mutations
//70    CCLE_ONCOMAP_total_mutations_in_gene
//71    CGC_Mutation_Type
//72    CGC_Translocation_Partner
//73    CGC_Tumor_Types_Somatic
//74    CGC_Tumor_Types_Germline
//75    CGC_Other_Diseases
//76    DNARepairGenes_Role
//77    FamilialCancerDatabase_Syndromes
//78    MUTSIG_Published_Results
//79    OREGANNO_ID
//80    OREGANNO_Values

}
