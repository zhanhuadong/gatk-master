import numpy as np
from gcnvkernel import Interval
from typing import Optional, Set, Generator
from intervaltree_bio import GenomeIntervalTree, IntervalTree


class GenericCopyNumberVariant(Interval):
    """This class represent a generic CNV locus.

    Note:
        Equality testing and hashing is done based similarly to Interval.
    """
    def __init__(self,
                 contig: str, start: int, end: int,
                 var_copy_number: int,
                 quality: float,
                 ref_copy_number: Optional[int] = None,
                 variant_frequency: Optional[float] = None,
                 num_intervals: Optional[int] = 1):
        """Initializer.

        Args:
            contig: contig
            start: start
            end: end
            var_copy_number: var copy number
            quality: quality score
            ref_copy_number: (optional) ref copy number
            variant_frequency: (optional) variant frequency in the cohort
            num_intervals: (optional) number of intervals that are merged to create the variant
        """
        super().__init__(contig, start, end)
        self.var_copy_number = var_copy_number
        self.ref_copy_number = ref_copy_number
        self.quality = quality
        self.variant_frequency = variant_frequency
        self.num_intervals = num_intervals

    @property
    def is_dup(self):
        assert self.ref_copy_number is not None, "Ref copy number metadata for {0} is not available. " \
                                                 "Cannot determine whether the variant is DUP.".format(self)
        return self.var_copy_number > self.ref_copy_number

    @property
    def is_del(self):
        assert self.ref_copy_number is not None, "Ref copy number metadata for {0} is not available. " \
                                                 "Cannot determine whether the variant is DEL.".format(self)
        return self.var_copy_number < self.ref_copy_number

    @property
    def is_var(self):
        assert self.ref_copy_number is not None, "Ref copy number metadata for {0} is not available. " \
                                                 "Cannot determine whether the variant is VAR.".format(self)
        return self.var_copy_number != self.ref_copy_number

    def variant_frequency_below_max_value(self, max_variant_frequency: float):
        assert self.variant_frequency is not None, "Variant {0} has no variant frequency annotation".format(self)
        return self.variant_frequency <= max_variant_frequency

    def quality_above_min_value(self, min_quality: float):
        return self.quality >= min_quality

    def num_intervals_above_min_value(self, min_num_intervals: int):
        assert self.variant_frequency is not None, "Variant {0} has no number of intervals annotation".format(self)
        return self.num_intervals >= min_num_intervals

    def get_padded(self, padding: int, keep_annotations=False):
        return GenericCopyNumberVariant(self.contig, self.start - padding, self.end + padding,
                                        self.var_copy_number, self.quality, self.ref_copy_number,
                                        self.variant_frequency, self.num_intervals)

    def get_symmetric_overlap_fraction(self, interval: Interval):
        if self.contig != interval.contig:
            return 0.0
        else:
            span = max(self.end, interval.end) - min(self.start, interval.start)
            span = span if span > 0 else 1e-12
            shared = max(0, min(self.end, interval.end) - max(self.start, interval.start))
            return float(shared) / span

    def __repr__(self):
        return "({0}, {1}, {2}), VAR_CN: {3}, REF_CN: {4}, QUAL: {5:.3f}, VF: {6}, NI: {7}".format(
            self.contig, self.start, self.end,
            self.var_copy_number,
            "N/A" if self.ref_copy_number is None else self.ref_copy_number,
            self.quality,
            "N/A" if self.variant_frequency is None else "{0:.3f}".format(self.variant_frequency),
            "N/A" if self.num_intervals is None else self.num_intervals)

    __str__ = __repr__


class GenericCNVCallSet:
    def __init__(self, sample_name: str):
        self.sample_name = sample_name
        self.genomic_interval_tree = GenomeIntervalTree()

    def add(self, variant: GenericCopyNumberVariant):
        self.genomic_interval_tree.addi(variant.contig, variant.start, variant.end, data=variant)

    def get_overlapping_variants_set(self,
                                     interval: Interval,
                                     min_symmetric_overlap_fraction: Optional[float] = None) \
            -> Set[GenericCopyNumberVariant]:
        tree_query = self.genomic_interval_tree[interval.contig].search(interval.start, interval.end)
        found_variants = {found_interval.data for found_interval in tree_query}
        if min_symmetric_overlap_fraction is None:
            return found_variants
        else:
            return {variant for variant in found_variants
                    if variant.get_symmetric_overlap_fraction(interval) >= min_symmetric_overlap_fraction}

    def overlaps_with(self, interval: Interval) -> bool:
        return self.genomic_interval_tree[interval.contig].overlaps(interval.start, interval.end)

    def iter_in_contig(self, contig: str) -> Generator:
        for entry in self.genomic_interval_tree[contig].iter():
            yield entry.data

    def get_contig_interval_tree(self, contig: str) -> IntervalTree:
        return self.genomic_interval_tree[contig]

    @property
    def contigs_set(self) -> Set[str]:
        return set(self.genomic_interval_tree.keys())

    def merge_nearby_variants(self, padding: int):
        """Returns a new call set in which nearby variants are merged together."""

        def get_unpadded_merged_variant(_overlapping_variants_set: Set[GenericCopyNumberVariant]) \
                -> GenericCopyNumberVariant:
            """From a set of overlapping variants, chooses the variant with highest quality as
            a representative, and returns a new variant that envelopes the overlapping set with
            the same attributes as the highest quality variant."""
            assert len(_overlapping_variants_set) > 0

            contig = next(iter(_overlapping_variants_set)).contig
            start = min([variant.start for variant in _overlapping_variants_set]) + padding
            end = max([variant.end for variant in _overlapping_variants_set]) - padding

            highest_quality_variant = sorted(_overlapping_variants_set, key=lambda variant: variant.quality)[-1]
            ref_copy_number = highest_quality_variant.ref_copy_number
            var_copy_number = highest_quality_variant.var_copy_number
            quality = highest_quality_variant.quality
            variant_frequency = highest_quality_variant.variant_frequency
            num_intervals = sum([variant.num_intervals for variant in _overlapping_variants_set])

            return GenericCopyNumberVariant(contig, start, end,
                                            var_copy_number, quality,
                                            ref_copy_number=ref_copy_number,
                                            variant_frequency=variant_frequency,
                                            num_intervals=num_intervals)

        # create an empty genome interval tree
        merged_genome_interval_tree = GenomeIntervalTree()

        # merge variants in each contig if necessary
        for contig in self.contigs_set:
            sorted_var_list = [iv.data for iv in sorted(list(self.get_contig_interval_tree(contig)))]

            merged_interval_tree = IntervalTree()
            current_overlapping_set = set()
            current_enveloping_interval = None
            for variant in sorted_var_list:
                padded_variant = variant.get_padded(padding)
                if len(current_overlapping_set) == 0:
                    current_overlapping_set.add(padded_variant)
                    current_enveloping_interval = Interval(
                        padded_variant.contig, padded_variant.start, padded_variant.end)
                elif current_enveloping_interval.overlaps_with(padded_variant):
                    current_overlapping_set.add(padded_variant)
                    current_enveloping_interval = Interval(
                        padded_variant.contig,
                        min(padded_variant.start, current_enveloping_interval.start),
                        max(padded_variant.end, current_enveloping_interval.end))
                else:
                    merged_variant = get_unpadded_merged_variant(current_overlapping_set)
                    merged_interval_tree.addi(merged_variant.start, merged_variant.end, data=merged_variant)
                    current_overlapping_set.clear()
                    current_overlapping_set.add(padded_variant)
                    current_enveloping_interval = Interval(
                        padded_variant.contig, padded_variant.start, padded_variant.end)
            if len(current_overlapping_set) > 0:
                merged_variant = get_unpadded_merged_variant(current_overlapping_set)
                merged_interval_tree.addi(merged_variant.start, merged_variant.end,
                                          data=merged_variant)

            merged_genome_interval_tree[contig] = merged_interval_tree

        merged_call_set = GenericCNVCallSet(self.sample_name)
        merged_call_set.genomic_interval_tree = merged_genome_interval_tree

        return merged_call_set



class TruthAndTrialVariants:
    """Stores a par of truth and trial variants."""
    def __init__(self,
                 truth_variant: GenericCopyNumberVariant,
                 trial_variant: GenericCopyNumberVariant):
        self.truth_variant = truth_variant
        self.trial_variant = trial_variant

    def __repr__(self):
        return "[truth: {0}\ntrial: {1}]".format(self.truth_variant, self.trial_variant)


class CNVCallSetAnalysisSummary:
    """Stores a summary of CNV call set analysis and provides a filtering mechanism."""

    def __init__(self):
        # overlaps with a truth variant, has the same copy number call
        self.exact_matches: Set[TruthAndTrialVariants] = set()

        # overlaps with a truth variant, copy number call is not the same as truth but both are dups
        self.qualitative_dup_matches: Set[TruthAndTrialVariants] = set()

        # overlaps with a truth variant, copy number call is not the same as truth but both are dels
        self.qualitative_del_matches: Set[TruthAndTrialVariants] = set()

        # truth is dup(del) but call is del(dup)
        self.opposite_calls: Set[TruthAndTrialVariants] = set()

        # dup truth variants that are not called (the set contains such truth variants)
        self.missed_dup_calls: Set[GenericCopyNumberVariant] = set()

        # del truth variants that are not called (the set contains such truth variants)
        self.missed_del_calls: Set[GenericCopyNumberVariant] = set()

        # truth is ref but call is dup (the set contains such trial variants)
        self.false_dup_calls: Set[GenericCopyNumberVariant] = set()

        # truth is ref but call is del (the set contains such trial variants)
        self.false_del_calls: Set[GenericCopyNumberVariant] = set()

        # trial variant excluded in truth analysis (the set contains such truth variants)
        self.trial_excluded_calls: Set[GenericCopyNumberVariant] = set()

        # truth variant excluded in trial analysis (the set contains such trial variants)
        self.truth_excluded_calls: Set[GenericCopyNumberVariant] = set()

    @property
    def true_positive_count(self):
        exact_matches_count = len(self.exact_matches)
        qualitative_dup_matches_count = len(self.qualitative_dup_matches)
        qualitative_del_matches_count = len(self.qualitative_del_matches)
        return exact_matches_count + qualitative_dup_matches_count + qualitative_del_matches_count

    @property
    def false_positive_count(self):
        false_dup_calls_count = len(self.false_dup_calls)
        false_del_calls_count = len(self.false_del_calls)
        return false_dup_calls_count + false_del_calls_count

    @property
    def false_negative_count(self):
        missed_dup_calls_count = len(self.missed_dup_calls)
        missed_del_calls_count = len(self.missed_del_calls)
        return missed_dup_calls_count + missed_del_calls_count

    @property
    def opposite_calls_count(self):
        return len(self.opposite_calls)

    @property
    def trial_excluded_calls_count(self):
        return len(self.trial_excluded_calls)

    @property
    def truth_excluded_calls_count(self):
        return len(self.truth_excluded_calls)

    @property
    def sensitivity(self):
        total_truth_calls = self.true_positive_count +  self.false_negative_count
        if total_truth_calls > 0:
            return float(self.true_positive_count) / total_truth_calls
        else:
            return np.nan

    @property
    def precision(self):
        total_trial_calls = self.true_positive_count + self.false_positive_count + self.opposite_calls_count
        if total_trial_calls > 0:
            return float(self.true_positive_count) / total_trial_calls
        else:
            return np.nan

    def get_filtered(self,
                     min_truth_quality: Optional[float] = None,
                     min_trial_quality: Optional[float] = None,
                     max_truth_variant_frequency: Optional[float] = None,
                     max_trial_variant_frequency: Optional[float] = None,
                     contig_set: Optional[Set[str]] = None):
        """Filters the summary according to the provided arguments."""

        def get_filtered_truth_and_trial_variants_set(input_set: Set[TruthAndTrialVariants]):
            output_set = set()
            for entry in input_set:
                if min_truth_quality is not None and entry.truth_variant.quality_above_min_value(min_truth_quality):
                    continue
                if min_trial_quality is not None and entry.trial_variant.quality_above_min_value(min_trial_quality):
                    continue
                if (max_truth_variant_frequency is not None and
                        entry.truth_variant.variant_frequency_below_max_value(max_truth_variant_frequency)):
                    continue
                if (max_trial_variant_frequency is not None and
                        entry.trial_variant.variant_frequency_below_max_value(max_trial_variant_frequency)):
                    continue
                if contig_set is not None and entry.truth_variant.contig not in contig_set:
                    continue
                output_set.add(entry)
            return output_set

        def get_filtered_truth_variants_set(input_set: Set[GenericCopyNumberVariant]):
            output_set = set()
            for entry in input_set:
                if min_truth_quality is not None and entry.quality_above_min_value(min_truth_quality):
                    continue
                if (max_truth_variant_frequency is not None and
                        entry.variant_frequency_below_max_value(max_truth_variant_frequency)):
                    continue
                if contig_set is not None and entry.contig not in contig_set:
                    continue
                output_set.add(entry)
            return output_set

        def get_filtered_trial_variants_set(input_set: Set[GenericCopyNumberVariant]):
            output_set = set()
            for entry in input_set:
                if min_trial_quality is not None and entry.quality_above_min_value(min_trial_quality):
                    continue
                if (max_trial_variant_frequency is not None and
                        entry.variant_frequency_below_max_value(max_trial_variant_frequency)):
                    continue
                if contig_set is not None and entry.contig not in contig_set:
                    continue
                output_set.add(entry)
            return output_set

        filtered_summary = CNVCallSetAnalysisSummary()

        filtered_summary.exact_matches = get_filtered_truth_and_trial_variants_set(
            self.exact_matches)

        filtered_summary.qualitative_dup_matches = get_filtered_truth_and_trial_variants_set(
            self.qualitative_dup_matches)

        filtered_summary.qualitative_del_matches = get_filtered_truth_and_trial_variants_set(
            self.qualitative_del_matches)

        filtered_summary.opposite_calls = get_filtered_truth_and_trial_variants_set(
            self.opposite_calls)

        filtered_summary.missed_dup_calls = get_filtered_truth_variants_set(
            self.missed_dup_calls)

        filtered_summary.missed_del_calls = get_filtered_truth_variants_set(
            self.missed_del_calls)

        filtered_summary.false_dup_calls = get_filtered_trial_variants_set(
            self.false_dup_calls)

        filtered_summary.false_del_calls = get_filtered_trial_variants_set(
            self.false_del_calls)

        filtered_summary.trial_excluded_calls = get_filtered_trial_variants_set(
            self.trial_excluded_calls)

        filtered_summary.truth_excluded_calls = get_filtered_truth_variants_set(
            self.truth_excluded_calls)

        return filtered_summary

    def print_summary_combined(self, prefix=""):
        print(prefix + "- Number of exact matches: {0}".format(len(self.exact_matches)))
        print(prefix + "- Number of qualitative DUP matches: {0}".format(len(self.qualitative_dup_matches)))
        print(prefix + "- Number of qualitative DEL matches: {0}".format(len(self.qualitative_del_matches)))
        print(prefix + "- Number of missed DUP variants: {0}".format(len(self.missed_dup_calls)))
        print(prefix + "- Number of missed DEL variants: {0}".format(len(self.missed_del_calls)))
        print(prefix + "- Number of false DUP variants: {0}".format(len(self.false_dup_calls)))
        print(prefix + "- Number of false DEL variants: {0}".format(len(self.false_del_calls)))
        print(prefix + "- Number of opposite calls: {0}".format(len(self.opposite_calls)))
        print(prefix + "- Number of excluded truth variants: {0}".format(len(self.truth_excluded_calls)))
        print(prefix + "- Number of excluded trial variants: {0}".format(len(self.trial_excluded_calls)))
        print(prefix + "* Sensitivity: {0:.3f}".format(self.sensitivity))
        print(prefix + "* Precision: {0:.3f}".format(self.precision))

    def print_summary_per_contig(self, contig_set: Set[str]):
        for contig in sorted(contig_set):
            print("Summary for contig {0}:".format(contig))
            print()
            self.get_filtered(contig_set={contig}).print_summary_combined("\t")
            print()


class CNVTrialCallSetEvaluator:
    def __init__(self,
                 truth_call_set: GenericCNVCallSet,
                 trial_call_set: GenericCNVCallSet,
                 truth_excluded_loci: GenomeIntervalTree = GenomeIntervalTree(),
                 trial_excluded_loci: GenomeIntervalTree = GenomeIntervalTree(),
                 min_overlap_fraction: float = 0.5,
                 exclusion_criterion: str = 'partial_overlap'):
        """Evaluates a trial CNV call set against a truth call set in genomic regions analyzed by both
        the truth and the trial callers.


        Args:
            truth_call_set: truth call set
            trial_call_set: trial call set
            truth_excluded_loci: genomic loci excluded by the truth caller
            trial_excluded_loci: genomic loci excluded by the trial caller
            min_overlap_fraction: consider two variants overlapping only if their symmetric overlap fraction
                is higher than this value
            exclusion_criterion: strategy to exclude a variant from either call set.
                "partial_overlap" results in a variant to be excluded even if it partially overlaps with
                an excluded interval. "full_overlap" requires full envelopment by excluded loci to be disqualified.
        """
        assert truth_call_set.sample_name == trial_call_set.sample_name,\
            "CallSets have different samples names; truth: {0}, trial: {1}".format(
                truth_call_set.sample_name, trial_call_set.sample_name)

        self.sample_name = truth_call_set.sample_name
        self.truth_call_set = truth_call_set
        self.trial_call_set = trial_call_set
        self.truth_excluded_loci = truth_excluded_loci
        self.trial_excluded_loci = trial_excluded_loci
        self.min_overlap_fraction = min_overlap_fraction
        self.all_contigs = truth_call_set.contigs_set.union(trial_call_set.contigs_set)

        if exclusion_criterion == "partial_overlap":
            self.strict = False
        elif exclusion_criterion == "full_overlap":
            self.strict = True
        else:
            raise Exception("Unknown exclusion criterion; choices: 'partial_overlap', 'full_overlap'")

    def __call__(self) -> CNVCallSetAnalysisSummary:
        """Generates the analysis summary."""

        def perform_analysis_for_contig(_contig: str, _summary: CNVCallSetAnalysisSummary) -> None:
            """Performs the evaluation for a given contig and updates the analysis summary."""
            truth_excluded_intervals = self.truth_excluded_loci[_contig]
            trial_excluded_intervals = self.trial_excluded_loci[_contig]

            # keeps track of truth variants that do not overlap with trial variants
            remaining_truth_variants = self.truth_call_set.get_contig_interval_tree(_contig).copy()

            # first, iterate over trial variants
            for trial_variant in self.trial_call_set.iter_in_contig(_contig):

                # if trial variant overlaps truth excluded intervals, neglect it
                if len(truth_excluded_intervals.search(
                        trial_variant.start, trial_variant.end, strict=self.strict)) > 0:
                    _summary.trial_excluded_calls.add(trial_variant)
                    continue

                # if trial variant overlaps with its own excluded intervals, neglect it
                if len(trial_excluded_intervals.search(
                        trial_variant.start, trial_variant.end, strict=self.strict)) > 0:
                    _summary.trial_excluded_calls.add(trial_variant)
                    continue

                # get a set of overlapping truth variants
                truth_overlapping_variants = self.truth_call_set.get_overlapping_variants_set(
                    trial_variant, self.min_overlap_fraction)

                if len(truth_overlapping_variants) == 0:  # not in truth

                    if trial_variant.is_dup:
                        _summary.false_dup_calls.add(trial_variant)
                    else:
                        _summary.false_del_calls.add(trial_variant)

                else:

                    # remove truth from the remaining truth variants
                    for truth_variant in truth_overlapping_variants:
                        try:
                            remaining_truth_variants.removei(truth_variant.start, truth_variant.end,
                                                             data=truth_variant)
                        except ValueError:  # a variant could be already removed
                            pass

                    # choose the highest quality truth variant from the overlapping set
                    best_truth_variant = max(truth_overlapping_variants, key=lambda variant: variant.quality)
                    joint_variant = TruthAndTrialVariants(best_truth_variant, trial_variant)

                    # exact match
                    if best_truth_variant.var_copy_number == trial_variant.var_copy_number:
                        _summary.exact_matches.add(joint_variant)
                        continue

                    # qualitative dup match
                    if best_truth_variant.is_dup and trial_variant.is_dup:
                        _summary.qualitative_dup_matches.add(joint_variant)
                        continue

                    # qualitative del match
                    if best_truth_variant.is_del and trial_variant.is_del:
                        _summary.qualitative_del_matches.add(joint_variant)
                        continue

                    # opposite
                    _summary.opposite_calls.add(joint_variant)

            # next, iterate over remaining truth variants
            for truth_entry in remaining_truth_variants.iter():
                truth_variant: GenericCopyNumberVariant = truth_entry.data

                # if truth variant is excluded in the trial analysis, neglect it
                if len(trial_excluded_intervals.search(
                        truth_variant.start, truth_variant.end, strict=self.strict)) > 0:
                    _summary.truth_excluded_calls.add(truth_variant)
                    continue

                # if truth variant is excluded by its own excluded loci, neglect it
                if len(truth_excluded_intervals.search(
                        truth_variant.start, truth_variant.end, strict=self.strict)) > 0:
                    _summary.truth_excluded_calls.add(truth_variant)
                    continue

                if truth_variant.is_dup:
                    _summary.missed_dup_calls.add(truth_variant)
                else:
                    _summary.missed_del_calls.add(truth_variant)

        summary = CNVCallSetAnalysisSummary()

        for contig in self.all_contigs:
            perform_analysis_for_contig(contig, summary)
        return summary
