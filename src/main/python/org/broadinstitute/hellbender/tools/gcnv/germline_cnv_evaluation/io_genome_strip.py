import vcf
from typing import List, Optional, Set
from .core import GenericCNVCallSet, GenericCopyNumberVariant
import logging

_logger = logging.getLogger(__name__)


def _get_mode_int_list(array: List[int]):
    """Finds the mode of an integer list."""
    _dict = {v: 0 for v in range(max(array) + 1)}
    for entry in array:
        _dict[entry] += 1
    return sorted([(k, v) for k, v in _dict.items()], key=lambda x: x[1])[-1][0]


def load_genome_strip_vcf_file(genome_strip_vcf_file: str,
                               max_records: Optional[int] = None,
                               log_frequency: int = 500,
                               allosomal_contigs: Set[str] = {'X', 'Y'},
                               autosomal_ref_copy_number: int = 2) -> List[GenericCNVCallSet]:
    """Loads Genome STRiP-style .VCF file and generates a list of call sets, one for each
    sample.

    Args:
        genome_strip_vcf_file: input Genome STRiP .VCF file
        max_records: (optional) maximum number of records to process
        log_frequency: log for every `log_frequency` processed records
        allosomal_contigs: set of allosomal contigs; used to determine ref copy number on such
            contigs
        autosomal_ref_copy_number: ref copy number for autosomal variants

    Note:
        The following replacements must be manually applied to any Genome STRiP VCF file in order
        to become compliant with VCF specs:
            - replace '=NA' with '=.' everywhere
            - replace '"NA"' with '"."' everywhere
            - replace '= "NA"' with '= "."' everywhere
            - change GSCNCATEGORY type from Float => String

    Returns:
        a list of call sets
    """

    # step 1. load VCF file
    with open(genome_strip_vcf_file, 'r') as f:
        vcf_reader = vcf.Reader(f)
        sample_names = vcf_reader.samples
        gs_call_set_list = list()
        for sample_name in sample_names:
            gs_call_set_list.append(GenericCNVCallSet(sample_name))

        for record_num, record in enumerate(vcf_reader):

            if record_num % log_frequency == 0:
                _logger.info("Reading record number {0}...".format(record_num))

            if max_records is not None and record_num > max_records:
                break

            info = record.INFO
            contig = record.CHROM
            start = record.start + 1
            end = info['END']

            for si, sample_record in enumerate(record.samples):
                # recalculate quality since the VCF caps it to 99
                # for some HOMDELs, CNP can be singleton; we set the quality to the segment quality
                copy_number_probs = sample_record.data.CNP
                quality = -10 * sorted(sample_record.data.CNP)[-2] \
                    if len(copy_number_probs) > 1 else info['GSCNQUAL']

                # generate variant
                sample_var = GenericCopyNumberVariant(contig, start, end, sample_record.data.CN,
                                                      quality, ref_copy_number=None, variant_frequency=None)
                gs_call_set_list[si].add(sample_var)

    # step 2. set the ref copy number
    #
    # Note:
    #   Since GS VCF contains variants discovered on any sample for all samples, the majority of
    #   variant contexts are REF for any given sample. We use this observation to determine the
    #   ploidy of sex chromosomes by identifying the most commonly called copy number state in
    #   sex chromosomes for each sample. For autosomal variants, we set `ref_copy_number` to
    #   `autosomal_ref_copy_number` (default=2).
    _logger.info("Determining reference copy numbers...")
    for gs_call_set in gs_call_set_list:
        _logger.info("Sample name: {0}".format(gs_call_set.sample_name))
        for contig in gs_call_set.contigs_set:
            if contig in allosomal_contigs:
                inferred_ref_for_contig = _get_mode_int_list(
                    [iv.data.var_copy_number for iv in gs_call_set.get_contig_interval_tree(contig)])
                _logger.info("contig: {0}, inferred REF_CN: {1}".format(contig, inferred_ref_for_contig))
            else:
                inferred_ref_for_contig = autosomal_ref_copy_number
            for iv in gs_call_set.get_contig_interval_tree(contig):
                iv.data.ref_copy_number = inferred_ref_for_contig

    # step 3. calculate variant frequencies
    _logger.info("Calculating variant frequencies...")
    contigs_set = gs_call_set_list[0].contigs_set
    for contig in contigs_set:
        generators = [gs_call_set.iter_in_contig(contig) for gs_call_set in gs_call_set_list]
        while True:
            try:
                variants = [generator.__next__() for generator in generators]
                assert len(set(variants)) == 1  # assert that the variants indeed come from the same interval
                num_var_samples = sum([variant.is_var for variant in variants])
                variant_frequency = float(num_var_samples) / len(variants)
                for variant in variants:
                    variant.variant_frequency = variant_frequency
            except StopIteration:
                break

    # step 4. non-variant "variants" are not needed anymore -- remove them from the interval tree
    _logger.info("Removing non-variant calls...")
    var_only_gs_call_set_list = list()
    for gs_call_set in gs_call_set_list:
        var_only_gs_call_set = GenericCNVCallSet(gs_call_set.sample_name)
        for contig in gs_call_set.contigs_set:
            for variant in gs_call_set.iter_in_contig(contig):
                if variant.is_var:
                    var_only_gs_call_set.add(variant)
        var_only_gs_call_set.tags.add("Removed non-variant calls")
        var_only_gs_call_set_list.append(var_only_gs_call_set)

    return var_only_gs_call_set_list
