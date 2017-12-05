import numpy as np
import pickle
from gcnvkernel import Interval
from typing import Optional, Set, Generator, List
from intervaltree_bio import GenomeIntervalTree, IntervalTree
from collections import namedtuple
import logging

_logger = logging.getLogger(__name__)


class GenericCopyNumberVariant(Interval):
    """This class represent a generic CNV locus.

    Note:
        Equality testing and hashing is done similarly to 'Interval', i.e. calls and other
        annotations are ignored.
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

    def variant_frequency_below_value(self, max_variant_frequency: float):
        assert self.variant_frequency is not None, "Variant {0} has no variant frequency annotation".format(self)
        return self.variant_frequency <= max_variant_frequency

    def variant_frequency_above_value(self, min_variant_frequency: float):
        assert self.variant_frequency is not None, "Variant {0} has no variant frequency annotation".format(self)
        return self.variant_frequency >= min_variant_frequency

    def quality_below_value(self, max_quality: float):
        return self.quality <= max_quality

    def quality_above_value(self, min_quality: float):
        return self.quality >= min_quality

    def length_below_value(self, max_length: int):
        return (self.end - self.start) <= max_length

    def length_above_value(self, min_length: int):
        return (self.end - self.start) >= min_length

    def get_padded(self, padding: int, keep_annotations=False):
        return GenericCopyNumberVariant(self.contig, self.start - padding, self.end + padding,
                                        self.var_copy_number, self.quality, self.ref_copy_number,
                                        self.variant_frequency, self.num_intervals)

    @staticmethod
    def get_variant_filter(min_quality: Optional[float] = None,
                           max_quality: Optional[float] = None,
                           min_variant_frequency: Optional[float] = None,
                           max_variant_frequency: Optional[float] = None,
                           min_length: Optional[int] = None,
                           max_length: Optional[int] = None):

        def min_quality_pass(variant: GenericCopyNumberVariant):
            return min_quality is None or variant.quality_above_value(min_quality)

        def max_quality_pass(variant: GenericCopyNumberVariant):
            return max_quality is None or variant.quality_below_value(max_quality)

        def min_variant_frequency_pass(variant: GenericCopyNumberVariant):
            return min_variant_frequency is None or variant.variant_frequency_above_value(min_variant_frequency)

        def max_variant_frequency_pass(variant: GenericCopyNumberVariant):
            return max_variant_frequency is None or variant.variant_frequency_below_value(max_variant_frequency)

        def min_length_pass(variant: GenericCopyNumberVariant):
            return min_length is None or variant.length_above_value(min_length)

        def max_length_pass(variant: GenericCopyNumberVariant):
            return max_length is None or variant.length_below_value(max_length)

        def all_pass(variant: GenericCopyNumberVariant):
            return (min_quality_pass(variant) and max_quality_pass(variant) and
                    min_variant_frequency_pass(variant) and max_variant_frequency_pass(variant) and
                    min_length_pass(variant) and max_length_pass(variant))

        return all_pass

    def get_symmetric_overlap_fraction(self, interval: Interval):
        if self.contig != interval.contig:
            return 0.0
        else:
            span = max(self.end, interval.end) - min(self.start, interval.start)
            span = span if span > 0 else 1e-12
            shared = max(0, min(self.end, interval.end) - max(self.start, interval.start))
            return float(shared) / span

    def __repr__(self):
        return "({0}, {1}, {2}), VAR_CN: {3}, REF_CN: {4}, GQ: {5:.3f}, VF: {6}, NI: {7}".format(
            self.contig, self.start, self.end,
            self.var_copy_number,
            "N/A" if self.ref_copy_number is None else self.ref_copy_number,
            self.quality,
            "N/A" if self.variant_frequency is None else "{0:.3f}".format(self.variant_frequency),
            "N/A" if self.num_intervals is None else self.num_intervals)

    __str__ = __repr__


GenericCNVCallSetPickleBundle = namedtuple(
    'GenericCNVCallSetPickleBundle', 'sample_name, tags, genome_interval_tree')


class GenericCNVCallSet:
    """This class represents a generic CNV call set and is used for standardizing the output of different
    tools for fair comparison."""
    def __init__(self, sample_name: str, tags: Set[str] = set()):
        self.sample_name = sample_name
        self.tags = tags
        self.genome_interval_tree = GenomeIntervalTree()

    @staticmethod
    def from_genome_interval_tree(sample_name: str,
                                  genome_interval_tree: GenomeIntervalTree,
                                  tags: Set[str] = set()):
        call_set = GenericCNVCallSet(sample_name, tags=tags)
        call_set.genome_interval_tree = genome_interval_tree
        return call_set

    @staticmethod
    def from_pickle(pickle_file: str) -> List['GenericCNVCallSet']:
        with open(pickle_file, 'rb') as f:
            unpickler = pickle.Unpickler(f)
            call_set_list = list()
            while True:
                try:
                    pickle_bundle = unpickler.load()
                    call_set_list.append(GenericCNVCallSet.from_genome_interval_tree(
                        pickle_bundle.sample_name, pickle_bundle.genome_interval_tree, pickle_bundle.tags))
                except EOFError:
                    break
            return call_set_list

    @staticmethod
    def to_pickle(pickle_file: str, call_set_list: List['GenericCNVCallSet']):
        with open(pickle_file, 'wb') as f:
            for call_set in call_set_list:
                pickle_bundle = GenericCNVCallSetPickleBundle(
                    call_set.sample_name,
                    call_set.tags,
                    call_set.genome_interval_tree)
                pickle.dump(pickle_bundle, f)

    def add(self, variant: GenericCopyNumberVariant):
        self.genome_interval_tree.addi(variant.contig, variant.start, variant.end, data=variant)

    def get_overlapping_variants_set(self,
                                     interval: Interval,
                                     min_symmetric_overlap_fraction: Optional[float] = None) \
            -> Set[GenericCopyNumberVariant]:
        tree_query = self.genome_interval_tree[interval.contig].search(interval.start, interval.end)
        found_variants = {found_interval.data for found_interval in tree_query}
        if min_symmetric_overlap_fraction is None:
            return found_variants
        else:
            return {variant for variant in found_variants
                    if variant.get_symmetric_overlap_fraction(interval) >= min_symmetric_overlap_fraction}

    def overlaps_with(self, interval: Interval) -> bool:
        return self.genome_interval_tree[interval.contig].overlaps(interval.start, interval.end)

    def iter_in_contig(self, contig: str) -> Generator:
        for entry in self.genome_interval_tree[contig].iter():
            yield entry.data

    def get_contig_interval_tree(self, contig: str) -> IntervalTree:
        return self.genome_interval_tree[contig]

    @property
    def contigs_set(self) -> Set[str]:
        return set(self.genome_interval_tree.keys())

    def merge_overlapping_variants(self, padding: int,
                                   interval_consensus_strategy: str = 'highest_quality',
                                   call_consensus_strategy: str = 'highest_quality') -> 'GenericCNVCallSet':
        """Returns a new call set in which nearby variants are merged together.

        Args:
            padding: (positive integer) the amount by which the intervals are to be padded
                symmetrically before detecting overlapping sets
            interval_consensus_strategy: strategy for obtaining the consensus interval.
                the choices are are follows:
                `highest_quality`: select the highest quality variant in the set
                `envelope`: yields the smallest interval that envelopes all variants in the
                overlapping set
            call_consensus_strategy: strategy for obtaining the consensus call (i.e.
                var copy number, quality, etc.). the choices are as follows:
                'highest_quality`: calls from the highest quality variant will be copied
                `average`: use average of the calls in the overlapping set
        """

        assert interval_consensus_strategy in {'highest_quality', 'envelope'}, \
            "Unrecognized interval consensus strategy"
        assert call_consensus_strategy in {'highest_quality', 'average'}, \
            "Unrecognized call consensus strategy"

        def generate_consensus_variant(_overlapping_variants_set: Set[GenericCopyNumberVariant]) \
                -> GenericCopyNumberVariant:
            """From a set of overlapping variants, chooses the variant with highest quality as
            a representative, and returns a new variant that envelopes the overlapping set with
            the same attributes as the highest quality variant.

            Args:
                _overlapping_variants_set: a set of overlapping `GenericCopyNumberVariant`, each
                    padded symmetrically by `padding`
            """
            assert len(_overlapping_variants_set) > 0

            _contig = next(iter(_overlapping_variants_set)).contig
            num_intervals = sum([_variant.num_intervals for _variant in _overlapping_variants_set])

            highest_quality_variant = sorted(_overlapping_variants_set, key=lambda _variant: _variant.quality)[-1]

            if interval_consensus_strategy == 'highest_quality':
                start = highest_quality_variant.start + padding
                end = highest_quality_variant.end - padding
            elif interval_consensus_strategy == 'envelope':
                start = min([_variant.start for _variant in _overlapping_variants_set]) + padding
                end = max([_variant.end for _variant in _overlapping_variants_set]) - padding
            else:
                raise Exception("Should not reach here")

            if call_consensus_strategy == 'highest_quality':
                ref_copy_number = highest_quality_variant.ref_copy_number
                var_copy_number = highest_quality_variant.var_copy_number
                quality = highest_quality_variant.quality
                variant_frequency = highest_quality_variant.variant_frequency
            elif call_consensus_strategy == 'average':  # todo
                raise NotImplementedError
            else:
                raise Exception("Should not reach here")

            return GenericCopyNumberVariant(
                _contig, start, end, var_copy_number, quality,
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
                    merged_variant = generate_consensus_variant(current_overlapping_set)
                    merged_interval_tree.addi(merged_variant.start, merged_variant.end, data=merged_variant)
                    current_overlapping_set.clear()
                    current_overlapping_set.add(padded_variant)
                    current_enveloping_interval = Interval(
                        padded_variant.contig, padded_variant.start, padded_variant.end)
            if len(current_overlapping_set) > 0:
                merged_variant = generate_consensus_variant(current_overlapping_set)
                merged_interval_tree.addi(merged_variant.start, merged_variant.end,
                                          data=merged_variant)

            merged_genome_interval_tree[contig] = merged_interval_tree

        tags = self.tags.copy()
        tags.add("Merged overlapping intervals with padding={0}, "
                 "interval_consensus_strategy={1}, "
                 "call_consensus_strategy={2}".format(
                    padding, interval_consensus_strategy, call_consensus_strategy))
        return GenericCNVCallSet.from_genome_interval_tree(
            self.sample_name, merged_genome_interval_tree, tags)

    def filter(self,
               min_quality: Optional[float] = None,
               max_quality: Optional[float] = None,
               min_variant_frequency: Optional[float] = None,
               max_variant_frequency: Optional[float] = None,
               min_length: Optional[int] = None,
               max_length: Optional[int] = None,
               included_contigs: Optional[Set[str]] = None,
               excluded_contigs: Optional[Set[str]] = None):
        final_tags = self.tags.copy()
        if included_contigs is not None:
            final_contigs_set = included_contigs.copy()
            final_tags.add("Included contigs: {0}".format(included_contigs))
        else:
            final_contigs_set = self.contigs_set.copy()
        if excluded_contigs is not None:
            final_contigs_set = final_contigs_set.difference(excluded_contigs)
            final_tags.add("Excluded contigs: {0}".format(excluded_contigs))
        final_contigs_set = final_contigs_set.intersection(self.contigs_set)

        if len(final_contigs_set) == 0:
            return GenericCNVCallSet(self.sample_name, final_tags)

        if min_quality is not None:
            final_tags.add("Minimum quality of kept variants: {0}".format(min_quality))

        if max_quality is not None:
            final_tags.add("Maximum quality of kept variants: {0}".format(max_quality))

        if min_variant_frequency is not None:
            final_tags.add("Minimum variant frequency of kept variants: {0}".format(min_variant_frequency))

        if max_variant_frequency is not None:
            final_tags.add("Maximum variant frequency of kept variants: {0}".format(max_variant_frequency))

        if min_length is not None:
            final_tags.add("Minimum length of kept variants: {0}".format(min_length))

        if max_length is not None:
            final_tags.add("Maximum length of kept variants: {0}".format(max_length))

        final_genome_interval_tree = GenomeIntervalTree()
        all_filters_pass = GenericCopyNumberVariant.get_variant_filter(
            min_quality=min_quality,
            max_quality=max_quality,
            min_variant_frequency=min_variant_frequency,
            max_variant_frequency=max_variant_frequency,
            min_length=min_length,
            max_length=max_length)

        from intervaltree_bio import Interval
        for contig in final_contigs_set:
            final_genome_interval_tree[contig] = IntervalTree(
                intervals=[Interval(variant.start, variant.end, variant)
                           for variant in filter(all_filters_pass, self.iter_in_contig(contig))])

        return GenericCNVCallSet.from_genome_interval_tree(
            self.sample_name, final_genome_interval_tree, final_tags)


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
    """Stores a summary of CNV call set analysis."""

    def __init__(self):
        # overlaps with a truth variant, has the same copy number call
        self.exact_matches: Set[TruthAndTrialVariants] = set()

        # overlaps with a truth variant, copy number call is not the same as truth but both are DUP
        self.qualitative_dup_matches: Set[TruthAndTrialVariants] = set()

        # overlaps with a truth variant, copy number call is not the same as truth but both are DEL
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
        total_truth_calls = self.true_positive_count + self.false_negative_count
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

    def get_summary_for_contig_set(self, contig_set: Set[str]) -> 'CNVCallSetAnalysisSummary':
        """Only keeps certain contigs."""

        def get_filtered_truth_and_trial_variants_set(input_set: Set[TruthAndTrialVariants]):
            return {pair for pair in filter(lambda _pair: _pair.truth_variant.contig in contig_set, input_set)}

        def get_filtered_truth_variants_set(input_set: Set[GenericCopyNumberVariant]):
            return {variant for variant in filter(lambda _variant: _variant.contig in contig_set, input_set)}

        def get_filtered_trial_variants_set(input_set: Set[GenericCopyNumberVariant]):
            return {variant for variant in filter(lambda _variant: _variant.contig in contig_set, input_set)}

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

    def print_summary_for_contig_set(self, contig_set: Set[str]):
        for contig in sorted(contig_set):
            print("Summary for contig {0}:".format(contig))
            print()
            self.get_summary_for_contig_set(contig_set={contig}).print_summary_combined("\t")
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
        if truth_call_set.sample_name != trial_call_set.sample_name:
            _logger.warning("CallSets have different samples names; truth: {0}, trial: {1}".format(
                truth_call_set.sample_name, trial_call_set.sample_name))

        self.sample_name = truth_call_set.sample_name
        self.truth_call_set = truth_call_set
        self.trial_call_set = trial_call_set
        self.truth_excluded_loci = truth_excluded_loci
        self.trial_excluded_loci = trial_excluded_loci
        self.min_overlap_fraction = min_overlap_fraction
        self.all_contigs = truth_call_set.contigs_set.union(trial_call_set.contigs_set)

        if exclusion_criterion == 'partial_overlap':
            self.strict = False
        elif exclusion_criterion == 'full_overlap':
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
