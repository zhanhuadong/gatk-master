package org.broadinstitute.hellbender.tools.funcotator.mafOutput;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.tools.funcotator.DataSourceFuncotationFactory;
import org.broadinstitute.hellbender.tools.funcotator.Funcotation;
import org.broadinstitute.hellbender.tools.funcotator.Funcotator;
import org.broadinstitute.hellbender.tools.funcotator.OutputRenderer;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;

/**
 * A Funcotator output renderer for writing to MAF files.
 * Complies with version 2.4 of the MAF standard found here:
 * https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4
 * Created by jonn on 12/5/17.
 */
public class MafOutputRenderer extends OutputRenderer {

    //==================================================================================================================
    // Public Static Members:

    //==================================================================================================================
    // Private Static Members:

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

    //==================================================================================================================
    // Constructors:

    public MafOutputRenderer(final Path outputFilePath,
                             final List<DataSourceFuncotationFactory> dataSources,
                             final LinkedHashMap<String, String> unaccountedForDefaultAnnotations,
                             final LinkedHashMap<String, String> unaccountedForOverrideAnnotations) {

        this.outputFilePath = outputFilePath;

        dataSourceFactories = dataSources;

        // Merge the annotations into our manualAnnotations:
        manualAnnotations = new LinkedHashMap<>();
        manualAnnotations.putAll(unaccountedForDefaultAnnotations);
        manualAnnotations.putAll(unaccountedForOverrideAnnotations);

        // Cache the manual annotation string so we can pass it easily into any Funcotations:
        manualAnnotationSerializedString = (manualAnnotations.size() != 0 ? String.join( FIELD_DELIMITER, manualAnnotations.values() ) + FIELD_DELIMITER : "");

        initializeDefaultMapWithKeys();
    }

    //==================================================================================================================
    // Override Methods:

    @Override
    public void open() {

    }

    @Override
    public void close() {

    }

    @Override
    public void write(final VariantContext variant, final List<Funcotation> funcotations) {

        //TODO: For this to work correctly you'll need to group funcotations by reference / alternate allele pairs.

        // First separate the funcotations into Gencode and Other:
        final List<GencodeFuncotation> gencodeFuncotations = new ArrayList<>(funcotations.size() / 2);
        final List<Funcotation> otherFuncotations = new ArrayList<>(funcotations.size() / 2);
        for ( final Funcotation funcotation : funcotations ) {
            if ( GencodeFuncotation.class.isAssignableFrom(funcotation.getClass()) ) {
                gencodeFuncotations.add((GencodeFuncotation) funcotation);
            }
            else {
                otherFuncotations.add(funcotation);
            }
        }

        // Create our output map:
        final LinkedHashMap<String, Object> outputMap = new LinkedHashMap<>(defaultMap);

        // Set values for unused fields:
        outputMap.put("Score", UNUSED_STRING);
        outputMap.put("BAM_File", UNUSED_STRING);

        // Go through Gencode funcotations as they're the basis for each line:
        for ( final GencodeFuncotation gencodeFuncotation : gencodeFuncotations ) {
            outputMap.put("Hugo_symbol", gencodeFuncotation.getHugoSymbol());
            outputMap.put("NCBI_Build", gencodeFuncotation.getNcbiBuild());
            outputMap.put("Chromosome", gencodeFuncotation.getChromosome());
            outputMap.put("Start_Position", gencodeFuncotation.getStart());
            outputMap.put("End_Position", gencodeFuncotation.getEnd());
            outputMap.put("Strand", gencodeFuncotation.getTranscriptStrand());
            outputMap.put("VariantClassification", gencodeFuncotation.getVariantClassification());
            outputMap.put("Variant_Type", gencodeFuncotation.getVariantType());
            outputMap.put("Reference_Allele", gencodeFuncotation.getRefAllele());
            outputMap.put("Tumor_Seq_Allele1", gencodeFuncotation.getTumorSeqAllele1());
            outputMap.put("Tumor_Seq_Allele2", gencodeFuncotation.getTumorSeqAllele2());
        }

        for ( final Funcotation funcotation : otherFuncotations ) {
            //TODO: For this to work you'll need to make annotations work as maps ala issue #3919
        }

    }


    //==================================================================================================================
    // Static Methods:

    //==================================================================================================================
    // Instance Methods:

    /**
     * Initializes {@link MafOutputRenderer#defaultMap} with the default keys for the columns in a MAF file.
     */
    void initializeDefaultMapWithKeys() {
        
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
        
        // TCGA required "optional" fields:
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
        defaultMap.put("COSMIC_overlapping_mutations",                UNKNOWN_VALUE_STRING);
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
        defaultMap.put("1000Genome_AA",                               UNKNOWN_VALUE_STRING);
        defaultMap.put("1000Genome_AC",                               UNKNOWN_VALUE_STRING);
        defaultMap.put("1000Genome_AF",                               UNKNOWN_VALUE_STRING);
        defaultMap.put("1000Genome_AFR_AF",                           UNKNOWN_VALUE_STRING);
        defaultMap.put("1000Genome_AMR_AF",                           UNKNOWN_VALUE_STRING);
        defaultMap.put("1000Genome_AN",                               UNKNOWN_VALUE_STRING);
        defaultMap.put("1000Genome_ASN_AF",                           UNKNOWN_VALUE_STRING);
        defaultMap.put("1000Genome_AVGPOST",                          UNKNOWN_VALUE_STRING);
        defaultMap.put("1000Genome_CIEND",                            UNKNOWN_VALUE_STRING);
        defaultMap.put("1000Genome_CIPOS",                            UNKNOWN_VALUE_STRING);
        defaultMap.put("1000Genome_END",                              UNKNOWN_VALUE_STRING);
        defaultMap.put("1000Genome_ERATE",                            UNKNOWN_VALUE_STRING);
        defaultMap.put("1000Genome_EUR_AF",                           UNKNOWN_VALUE_STRING);
        defaultMap.put("1000Genome_HOMLEN",                           UNKNOWN_VALUE_STRING);
        defaultMap.put("1000Genome_HOMSEQ",                           UNKNOWN_VALUE_STRING);
        defaultMap.put("1000Genome_LDAF",                             UNKNOWN_VALUE_STRING);
        defaultMap.put("1000Genome_RSQ",                              UNKNOWN_VALUE_STRING);
        defaultMap.put("1000Genome_SNPSOURCE",                        UNKNOWN_VALUE_STRING);
        defaultMap.put("1000Genome_SVLEN",                            UNKNOWN_VALUE_STRING);
        defaultMap.put("1000Genome_SVTYPE",                           UNKNOWN_VALUE_STRING);
        defaultMap.put("1000Genome_THETA",                            UNKNOWN_VALUE_STRING);
        defaultMap.put("1000Genome_VT",                               UNKNOWN_VALUE_STRING);
        defaultMap.put("ACHILLES_Lineage_Results_Top_Genes",          UNKNOWN_VALUE_STRING);
        defaultMap.put("CGC_Cancer Germline Mut",                     UNKNOWN_VALUE_STRING);
        defaultMap.put("CGC_Cancer Molecular Genetics",               UNKNOWN_VALUE_STRING);
        defaultMap.put("CGC_Cancer Somatic Mut",                      UNKNOWN_VALUE_STRING);
        defaultMap.put("CGC_Cancer Syndrome",                         UNKNOWN_VALUE_STRING);
        defaultMap.put("CGC_Chr",                                     UNKNOWN_VALUE_STRING);
        defaultMap.put("CGC_Chr Band",                                UNKNOWN_VALUE_STRING);
        defaultMap.put("CGC_GeneID",                                  UNKNOWN_VALUE_STRING);
        defaultMap.put("CGC_Name",                                    UNKNOWN_VALUE_STRING);
        defaultMap.put("CGC_Other Germline Mut",                      UNKNOWN_VALUE_STRING);
        defaultMap.put("CGC_Tissue Type",                             UNKNOWN_VALUE_STRING);
        defaultMap.put("COSMIC_n_overlapping_mutations",              UNKNOWN_VALUE_STRING);
        defaultMap.put("COSMIC_overlapping_mutation_descriptions",    UNKNOWN_VALUE_STRING);
        defaultMap.put("COSMIC_overlapping_primary_sites",            UNKNOWN_VALUE_STRING);
        defaultMap.put("ClinVar_ASSEMBLY",                            UNKNOWN_VALUE_STRING);
        defaultMap.put("ClinVar_HGMD_ID",                             UNKNOWN_VALUE_STRING);
        defaultMap.put("ClinVar_SYM",                                 UNKNOWN_VALUE_STRING);
        defaultMap.put("ClinVar_TYPE",                                UNKNOWN_VALUE_STRING);
        defaultMap.put("ClinVar_rs",                                  UNKNOWN_VALUE_STRING);
        defaultMap.put("ESP_AA",                                      UNKNOWN_VALUE_STRING);
        defaultMap.put("ESP_AAC",                                     UNKNOWN_VALUE_STRING);
        defaultMap.put("ESP_AA_AC",                                   UNKNOWN_VALUE_STRING);
        defaultMap.put("ESP_AA_AGE",                                  UNKNOWN_VALUE_STRING);
        defaultMap.put("ESP_AA_GTC",                                  UNKNOWN_VALUE_STRING);
        defaultMap.put("ESP_AvgAAsampleReadDepth",                    UNKNOWN_VALUE_STRING);
        defaultMap.put("ESP_AvgEAsampleReadDepth",                    UNKNOWN_VALUE_STRING);
        defaultMap.put("ESP_AvgSampleReadDepth",                      UNKNOWN_VALUE_STRING);
        defaultMap.put("ESP_CA",                                      UNKNOWN_VALUE_STRING);
        defaultMap.put("ESP_CDP",                                     UNKNOWN_VALUE_STRING);
        defaultMap.put("ESP_CG",                                      UNKNOWN_VALUE_STRING);
        defaultMap.put("ESP_CP",                                      UNKNOWN_VALUE_STRING);
        defaultMap.put("ESP_Chromosome",                              UNKNOWN_VALUE_STRING);
        defaultMap.put("ESP_DBSNP",                                   UNKNOWN_VALUE_STRING);
        defaultMap.put("ESP_DP",                                      UNKNOWN_VALUE_STRING);
        defaultMap.put("ESP_EA_AC",                                   UNKNOWN_VALUE_STRING);
        defaultMap.put("ESP_EA_AGE",                                  UNKNOWN_VALUE_STRING);
        defaultMap.put("ESP_EA_GTC",                                  UNKNOWN_VALUE_STRING);
        defaultMap.put("ESP_EXOME_CHIP",                              UNKNOWN_VALUE_STRING);
        defaultMap.put("ESP_FG",                                      UNKNOWN_VALUE_STRING);
        defaultMap.put("ESP_GL",                                      UNKNOWN_VALUE_STRING);
        defaultMap.put("ESP_GM",                                      UNKNOWN_VALUE_STRING);
        defaultMap.put("ESP_GS",                                      UNKNOWN_VALUE_STRING);
        defaultMap.put("ESP_GTC",                                     UNKNOWN_VALUE_STRING);
        defaultMap.put("ESP_GTS",                                     UNKNOWN_VALUE_STRING);
        defaultMap.put("ESP_GWAS_PUBMED",                             UNKNOWN_VALUE_STRING);
        defaultMap.put("ESP_MAF",                                     UNKNOWN_VALUE_STRING);
        defaultMap.put("ESP_PH",                                      UNKNOWN_VALUE_STRING);
        defaultMap.put("ESP_PP",                                      UNKNOWN_VALUE_STRING);
        defaultMap.put("ESP_Position",                                UNKNOWN_VALUE_STRING);
        defaultMap.put("ESP_TAC",                                     UNKNOWN_VALUE_STRING);
        defaultMap.put("ESP_TotalAAsamplesCovered",                   UNKNOWN_VALUE_STRING);
        defaultMap.put("ESP_TotalEAsamplesCovered",                   UNKNOWN_VALUE_STRING);
        defaultMap.put("ESP_TotalSamplesCovered",                     UNKNOWN_VALUE_STRING);
        defaultMap.put("Ensembl_so_accession",                        UNKNOWN_VALUE_STRING);
        defaultMap.put("Ensembl_so_term",                             UNKNOWN_VALUE_STRING);
        defaultMap.put("Familial_Cancer_Genes_Reference",             UNKNOWN_VALUE_STRING);
        defaultMap.put("Familial_Cancer_Genes_Synonym",               UNKNOWN_VALUE_STRING);
        defaultMap.put("HGNC_Ensembl Gene ID",                        UNKNOWN_VALUE_STRING);
        defaultMap.put("HGNC_HGNC ID",                                UNKNOWN_VALUE_STRING);
        defaultMap.put("HGNC_RefSeq IDs",                             UNKNOWN_VALUE_STRING);
        defaultMap.put("HGNC_Status",                                 UNKNOWN_VALUE_STRING);
        defaultMap.put("HGNC_UCSC ID(supplied by UCSC)",              UNKNOWN_VALUE_STRING);
        defaultMap.put("HGVS_coding_DNA_change",                      UNKNOWN_VALUE_STRING);
        defaultMap.put("HGVS_genomic_change",                         UNKNOWN_VALUE_STRING);
        defaultMap.put("HGVS_protein_change",                         UNKNOWN_VALUE_STRING);
        defaultMap.put("ORegAnno_bin",                                UNKNOWN_VALUE_STRING);
        defaultMap.put("UniProt_alt_uniprot_accessions",              UNKNOWN_VALUE_STRING);
        defaultMap.put("build",                                       UNKNOWN_VALUE_STRING);
        defaultMap.put("ccds_id",                                     UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_1000Gp1_AC",                           UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_1000Gp1_AF",                           UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_1000Gp1_AFR_AC",                       UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_1000Gp1_AFR_AF",                       UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_1000Gp1_AMR_AC",                       UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_1000Gp1_AMR_AF",                       UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_1000Gp1_ASN_AC",                       UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_1000Gp1_ASN_AF",                       UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_1000Gp1_EUR_AC",                       UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_1000Gp1_EUR_AF",                       UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_Ancestral_allele",                     UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_CADD_phred",                           UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_CADD_raw",                             UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_CADD_raw_rankscore",                   UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_ESP6500_AA_AF",                        UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_ESP6500_EA_AF",                        UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_Ensembl_geneid",                       UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_Ensembl_transcriptid",                 UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_FATHMM_pred",                          UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_FATHMM_rankscore",                     UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_FATHMM_score",                         UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_GERP++_NR",                            UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_GERP++_RS",                            UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_GERP++_RS_rankscore",                  UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_Interpro_domain",                      UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_LRT_Omega",                            UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_LRT_converted_rankscore",              UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_LRT_pred",                             UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_LRT_score",                            UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_LR_pred",                              UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_LR_rankscore",                         UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_LR_score",                             UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_MutationAssessor_pred",                UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_MutationAssessor_rankscore",           UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_MutationAssessor_score",               UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_MutationTaster_converted_rankscore",   UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_MutationTaster_pred",                  UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_MutationTaster_score",                 UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_Polyphen2_HDIV_pred",                  UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_Polyphen2_HDIV_rankscore",             UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_Polyphen2_HDIV_score",                 UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_Polyphen2_HVAR_pred",                  UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_Polyphen2_HVAR_rankscore",             UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_Polyphen2_HVAR_score",                 UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_RadialSVM_pred",                       UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_RadialSVM_rankscore",                  UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_RadialSVM_score",                      UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_Reliability_index",                    UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_SIFT_converted_rankscore",             UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_SIFT_pred",                            UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_SIFT_score",                           UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_SLR_test_statistic",                   UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_SiPhy_29way_logOdds",                  UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_SiPhy_29way_logOdds_rankscore",        UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_SiPhy_29way_pi",                       UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_UniSNP_ids",                           UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_Uniprot_aapos",                        UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_Uniprot_acc",                          UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_Uniprot_id",                           UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_aaalt",                                UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_aapos",                                UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_aapos_FATHMM",                         UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_aapos_SIFT",                           UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_aaref",                                UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_cds_strand",                           UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_codonpos",                             UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_fold-degenerate",                      UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_genename",                             UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_hg18_pos(1-coor)",                     UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_phastCons100way_vertebrate",           UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_phastCons100way_vertebrate_rankscore", UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_phastCons46way_placental",             UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_phastCons46way_placental_rankscore",   UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_phastCons46way_primate",               UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_phastCons46way_primate_rankscore",     UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_phyloP100way_vertebrate",              UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_phyloP100way_vertebrate_rankscore",    UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_phyloP46way_placental",                UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_phyloP46way_placental_rankscore",      UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_phyloP46way_primate",                  UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_phyloP46way_primate_rankscore",        UNKNOWN_VALUE_STRING);
        defaultMap.put("dbNSFP_refcodon",                             UNKNOWN_VALUE_STRING);
        defaultMap.put("gencode_transcript_name",                     UNKNOWN_VALUE_STRING);
        defaultMap.put("gencode_transcript_status",                   UNKNOWN_VALUE_STRING);
        defaultMap.put("gencode_transcript_tags",                     UNKNOWN_VALUE_STRING);
        defaultMap.put("gencode_transcript_type",                     UNKNOWN_VALUE_STRING);
        defaultMap.put("gene_id",                                     UNKNOWN_VALUE_STRING);
        defaultMap.put("gene_type",                                   UNKNOWN_VALUE_STRING);
        defaultMap.put("havana_transcript",                           UNKNOWN_VALUE_STRING);
        defaultMap.put("secondary_variant_classification",            UNKNOWN_VALUE_STRING);
        defaultMap.put("strand",                                      UNKNOWN_VALUE_STRING);
        defaultMap.put("transcript_id",                               UNKNOWN_VALUE_STRING);
    }

    //==================================================================================================================
    // Helper Data Types:

    //------------------------------------------------------------------------------------------------------------------
    // Required Columns:

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

}
