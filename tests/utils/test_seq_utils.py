from unittest import TestCase

from make_prg import MSA
from make_prg.utils.misc import equal_msas
from make_prg.utils.seq_utils import (
    GAP,
    NONMATCH,
    SequenceCurationError,
    SequenceExpander,
    align,
    count,
    get_alignment_seqs,
    get_consensus_from_MSA,
    get_majority_consensus_from_MSA,
    get_number_of_unique_gapped_sequences,
    get_number_of_unique_ungapped_sequences,
    has_empty_sequence,
    is_gap,
    is_non_match,
    remove_columns_full_of_gaps_from_MSA,
    remove_duplicates,
    ungap,
)
from tests.test_helpers import make_alignment


class TestSeqUtilsMisc(TestCase):
    def test___is_non_match___column_is_all_As(self):
        expected = False
        actual = is_non_match("A")
        self.assertEqual(expected, actual)

    def test___is_non_match___column_is_non_match(self):
        expected = True
        actual = is_non_match(NONMATCH)
        self.assertEqual(expected, actual)

    def test___is_gap___column_is_all_As(self):
        expected = False
        actual = is_gap("A")
        self.assertEqual(expected, actual)

    def test___is_gap___column_is_gap(self):
        expected = True
        actual = is_gap(GAP)
        self.assertEqual(expected, actual)

    def test___ungap___seq_has_no_gaps(self):
        sequence = "ACGTGTGACA"
        expected = "ACGTGTGACA"
        actual = ungap(sequence)
        self.assertEqual(expected, actual)

    def test___ungap___seq_has_gaps(self):
        sequence = "---AC-GT---GTG--AC--A"
        expected = "ACGTGTGACA"
        actual = ungap(sequence)
        self.assertEqual(expected, actual)

    def test___get_number_of_unique_ungapped_sequences___single_spaces(self):
        msa = make_alignment(
            [
                "-ACGT",
                "A-CGT",
                "AC-GT",
                "ACG-T",
                "ACGT-",
                "AC--T",
                "A-CGT",
                "A-CGT",
                "ACG-T",
                "ACG-T",
                "AC--T",
                "AC--T",
                "AC--T",
                "A-CGT",
            ]
        )

        expected = 2
        actual = get_number_of_unique_ungapped_sequences(msa)

        self.assertEqual(expected, actual)

    def test___get_number_of_unique_ungapped_sequences___several_spaces(self):
        msa = make_alignment(
            [
                "TA--------T",
                "T-A-------T",
                "T--A------T",
                "T---A-----T",
                "T----A----T",
                "T-----A---T",
                "T------A--T",
                "T-------A-T",
                "T--------AT",
                "-T-------AT",
                "--T------AT",
                "---T-----AT",
                "T-------AT-",
                "T------AT--",
                "T-----AT---",
                "T--A------T",
                "T--A------T",
                "---T-----AT",
                "---T-----AT",
                "--T------AT",
                "-T-------AT",
                "T-----A---T",
                "T------A--T",
                "T-------A-T",
            ]
        )

        expected = 1
        actual = get_number_of_unique_ungapped_sequences(msa)

        self.assertEqual(expected, actual)

    def test___get_number_of_unique_gapped_sequences___single_spaces(self):
        msa = make_alignment(
            [
                "-ACGT",
                "A-CGT",
                "AC-GT",
                "ACG-T",
                "ACGT-",
                "AC--T",
                "A-CGT",
                "A-CGT",
                "ACG-T",
                "ACG-T",
                "AC--T",
                "AC--T",
                "AC--T",
                "A-CGT",
            ]
        )

        expected = 6
        actual = get_number_of_unique_gapped_sequences(msa)

        self.assertEqual(expected, actual)

    def test___get_number_of_unique_gapped_sequences___several_spaces(self):
        msa = make_alignment(
            [
                "TA--------T",
                "T-A-------T",
                "T--A------T",
                "T---A-----T",
                "T----A----T",
                "T-----A---T",
                "T------A--T",
                "T-------A-T",
                "T--------AT",
                "-T-------AT",
                "--T------AT",
                "---T-----AT",
                "T-------AT-",
                "T------AT--",
                "T-----AT---",
                "T--A------T",
                "T--A------T",
                "---T-----AT",
                "---T-----AT",
                "--T------AT",
                "-T-------AT",
                "T-----A---T",
                "T------A--T",
                "T-------A-T",
            ]
        )

        expected = 15
        actual = get_number_of_unique_gapped_sequences(msa)

        self.assertEqual(expected, actual)


class TestRemoveDuplicates(TestCase):
    def test___remove_duplicates___empty_list(self):
        the_list = []

        expected = []
        actual = list(remove_duplicates(the_list))

        self.assertEqual(expected, actual)

    def test___remove_duplicates___unique_list(self):
        the_list = ["1", "2", "3"]

        expected = ["1", "2", "3"]
        actual = list(remove_duplicates(the_list))

        self.assertEqual(expected, actual)

    def test___remove_duplicates___several_copies(self):
        the_list = ["1", "1", "2", "3", "2", "2", "2", "2", "1", "1"]

        expected = ["1", "2", "3"]
        actual = list(remove_duplicates(the_list))

        self.assertEqual(expected, actual)


class TestEmptySeq(TestCase):
    def test___has_empty_sequence___individual_columns(self):
        msa = make_alignment(["ACCCT", "A---T", "A-G-T"])
        self.assertFalse(has_empty_sequence(msa, (0, 0)))
        self.assertTrue(has_empty_sequence(msa, (1, 1)))
        self.assertTrue(has_empty_sequence(msa, (2, 2)))
        self.assertTrue(has_empty_sequence(msa, (3, 3)))
        self.assertFalse(has_empty_sequence(msa, (4, 4)))

    def test___has_empty_sequence___long_columns(self):
        msa = make_alignment(["AAACCCGGGTTT", "AAA-C----T--", "AAA----G---T"])
        self.assertFalse(has_empty_sequence(msa, (0, 2)))
        self.assertTrue(has_empty_sequence(msa, (3, 5)))
        self.assertTrue(has_empty_sequence(msa, (6, 8)))
        self.assertFalse(has_empty_sequence(msa, (9, 11)))


class TestSequenceExpander(TestCase):
    def test___get_expanded_sequences_from_MSA___ambiguous_bases_one_seq(self):
        msa = make_alignment(["RWAAT"])

        expected = {"GAAAT", "AAAAT", "GTAAT", "ATAAT"}
        actual = SequenceExpander.get_expanded_sequences_from_MSA(msa)

        self.assertEqual(expected, set(actual))

    def test___get_expanded_sequences_from_MSA___ambiguous_bases_one_seq_with_repeated_base(
        self,
    ):
        msa = make_alignment(["RRAAT"])

        expected = {"GAAAT", "AAAAT", "GGAAT", "AGAAT"}
        actual = SequenceExpander.get_expanded_sequences_from_MSA(msa)

        self.assertEqual(expected, set(actual))

    def test___get_expanded_sequences_from_MSA___ambiguous_bases_several_seqs_with_gaps___all_wildcards(
        self,
    ):
        msa = make_alignment(
            [
                "--RYAA",
                "KMTT--",
                "-CSWG-",
                "GTAA--",  # one of the expanded ones
                "--ATAA",  # one of the expanded ones
                "-GCTT-",  # one of the expanded ones
                "CG--TG",  # one of the expanded ones
                "C-CT-G",  # one of the expanded ones
                "AA--AA",  # new non-wildcard seq
            ]
        )

        expected = {
            "GTAA",
            "GCAA",
            "ATAA",
            "ACAA",
            "GATT",
            "GCTT",
            "TATT",
            "TCTT",
            "CGAG",
            "CGTG",
            "CCAG",
            "CCTG",
            "AAAA",
        }
        actual = SequenceExpander.get_expanded_sequences_from_MSA(msa)

        self.assertEqual(expected, set(actual))

    def test___get_expanded_sequences_from_MSA___ambiguous_bases_several_seqs_with_gaps_and_Ns___all_wildcards(
        self,
    ):
        msa = make_alignment(
            ["-NRYAA", "KMTT--", "-CSWGN", "AA--AA"]  # new non-wildcard seq
        )

        expected = {"GATT", "GCTT", "TATT", "TCTT", "AAAA"}
        actual = SequenceExpander.get_expanded_sequences_from_MSA(msa)

        self.assertEqual(expected, set(actual))

    def test___get_expanded_sequences_from_MSA___ambiguous_bases_all_seqs_with_Ns___raises_SequenceCurationError(
        self,
    ):
        msa = make_alignment(
            ["-NRYAA", "KMNT--", "-CSWGN", "AA--NA"]  # new non-wildcard seq
        )

        with self.assertRaises(SequenceCurationError):
            SequenceExpander.get_expanded_sequences_from_MSA(msa)

    def test___get_expanded_sequences_from_MSA___first_sequence_in_is_first_sequence_out(
        self,
    ):
        msa = make_alignment(["TTTT", "AAAA", "CC-C"])

        expected = ["TTTT", "AAAA", "CCC"]
        actual = SequenceExpander.get_expanded_sequences_from_MSA(msa)

        self.assertEqual(expected, actual)

    def test___get_expanded_sequences_from_MSA___sequence_with_disallowed_char(self):
        msa = make_alignment(["TTTTN", "AAAAN", "CC-C-", "GGGGN", "AAAAB"])
        with self.assertRaises(SequenceCurationError):
            SequenceExpander.get_expanded_sequences_from_MSA(msa)


class TestAlign(TestCase):
    def test___align___small_sequence_with_indel___prioritise_matches(self):
        seq1 = "TA"
        seq2 = "TGA"

        expected = ("T-A", "TGA")
        actual = align(seq1, seq2)

        self.assertEqual(expected, actual)

    def test___align___medium_sequence_with_indel___prioritise_matches(self):
        seq1 = "TTTAAA"
        seq2 = "TTTGGAAA"

        expected = ("TTT--AAA", "TTTGGAAA")
        actual = align(seq1, seq2)

        self.assertEqual(expected, actual)

    def test___align___small_sequence_with_indel_and_mismatch(self):
        seq1 = "CA"
        seq2 = "TGA"

        expected = ("C-A", "TGA")
        actual = align(seq1, seq2)

        self.assertEqual(expected, actual)

    def test___align___several_matches_islands_with_indels___prioritise_matches(self):
        seq1 = "AAACCCGGGTTT"
        seq2 = "GTGAAAGGCCCTATAGGGAAATTTAA"

        expected = ("---AAA--CCC----GGG---TTT--", "GTGAAAGGCCCTATAGGGAAATTTAA")
        actual = align(seq1, seq2)

        self.assertEqual(expected, actual)

    def test___align___several_matches_islands_with_indels_and_mismatches___prioritise_matches(
        self,
    ):
        seq1 = "AAAAACCCCCGGGGGTTTTT"
        seq2 = "GTGAATAAGGCCGCCTATAGGCGGAAATTATTAA"

        expected = (
            "---AAAAA--CCCCC----GGGGG---TTTTT--",
            "GTGAATAAGGCCGCCTATAGGCGGAAATTATTAA",
        )
        actual = align(seq1, seq2)

        self.assertEqual(expected, actual)

    def test___align___seq1_is_empty(self):
        seq1 = ""
        seq2 = "ACGT"

        expected = ("----", "ACGT")
        actual = align(seq1, seq2)

        self.assertEqual(expected, actual)

    def test___align___seq2_is_empty(self):
        seq1 = "ACGT"
        seq2 = ""

        expected = ("ACGT", "----")
        actual = align(seq1, seq2)

        self.assertEqual(expected, actual)

    def test___align___both_seqs_are_empty(self):
        expected = ("", "")
        actual = align("", "")
        self.assertEqual(expected, actual)

    def test___align___several_equally_good_alignments___only_one_is_chosen(self):
        seq1 = "A"
        seq2 = "T"

        expected = ("A-", "-T")
        actual = align(seq1, seq2, 0, 0, 0, 0)

        self.assertEqual(expected, actual)


class Test_remove_columns_full_of_gaps_from_MSA(TestCase):
    def test___remove_columns_full_of_gaps_from_MSA___no_gaps_at_all(self):
        msa = make_alignment(
            [
                "AAAA",
                "AACA",
                "AAGA",
                "AATA",
            ]
        )

        expected = msa
        actual = remove_columns_full_of_gaps_from_MSA(msa)

        self.assertTrue(equal_msas(expected, actual))

    def test___remove_columns_full_of_gaps_from_MSA___some_gaps___no_columns_removed(
        self,
    ):
        msa = make_alignment(
            [
                "AAA-",
                "-ACA",
                "-AG-",
                "A-T-",
            ]
        )

        expected = msa
        actual = remove_columns_full_of_gaps_from_MSA(msa)

        self.assertTrue(equal_msas(expected, actual))

    def test___remove_columns_full_of_gaps_from_MSA___some_gaps___one_column_removed(
        self,
    ):
        msa = make_alignment(
            [
                "AA--",
                "-A-A",
                "-A--",
                "A---",
            ]
        )

        expected = make_alignment(
            [
                "AA-",
                "-AA",
                "-A-",
                "A--",
            ]
        )
        actual = remove_columns_full_of_gaps_from_MSA(msa)

        self.assertTrue(equal_msas(expected, actual))

    def test___remove_columns_full_of_gaps_from_MSA___some_gaps___two_columns_removed(
        self,
    ):
        msa = make_alignment(
            [
                "-A--",
                "-A-A",
                "-A--",
                "----",
            ]
        )

        expected = make_alignment(
            [
                "A-",
                "AA",
                "A-",
                "--",
            ]
        )
        actual = remove_columns_full_of_gaps_from_MSA(msa)

        self.assertTrue(equal_msas(expected, actual))


class TestSeqIteration(TestCase):
    input_seqs = ["TTAGGTTT", "TTA--TTT", "GGA-TTTT"]

    def test___get_seqs_from_alignment(self):
        msa = make_alignment(self.input_seqs)

        expected = self.input_seqs
        actual = list(get_alignment_seqs(msa))

        self.assertEqual(expected, actual)

    def test___count_alignment_seqs(self):
        msa = make_alignment(self.input_seqs)

        expected = 3
        actual = count(get_alignment_seqs(msa))

        self.assertEqual(expected, actual)


class Test_get_consensus_from_MSA(TestCase):
    def test___get_consensus_from_MSA___all_match(self):
        alignment = make_alignment(["AATTA", "AATTA"])

        expected = "AATTA"
        actual = get_consensus_from_MSA(alignment)

        self.assertEqual(expected, actual)

    def test___get_consensus_from_MSA___mixed_match_nonmatch(self):
        alignment = make_alignment(["AAGTA", "CATTA"])

        expected = "*A*TA"
        actual = get_consensus_from_MSA(alignment)

        self.assertEqual(expected, actual)

    def test___get_consensus_from_MSA___indel_nonmatch(self):
        alignment = make_alignment(["AAAA", "A--A"])

        expected = "A**A"
        actual = get_consensus_from_MSA(alignment)

        self.assertEqual(expected, actual)

    def test___get_consensus_from_MSA___IUPACAmbiguous_nonmatch(self):
        alignment = make_alignment(["RYA", "RTA"])

        expected = "**A"
        actual = get_consensus_from_MSA(alignment)

        self.assertEqual(expected, actual)

    def test___get_consensus_from_MSA___N_special_treatment(self):
        """
        i)A and N at pos 2 are different, but still consensus
        ii)N and N at pos 0 are same, but not consensus"""
        alignment = make_alignment(["TNNNN", "TNAR-"])

        expected = "T*A**"
        actual = get_consensus_from_MSA(alignment)

        self.assertEqual(expected, actual)

    def test___get_consensus_from_MSA___N_special_treatment___several_sequences(self):
        """
        i)A and N at pos 2 are different, but still consensus
        ii)N and N at pos 0 are same, but not consensus"""
        alignment = make_alignment(["TNNNNN", "TANAAN", "TAA-RN"])

        expected = "TAA***"
        actual = get_consensus_from_MSA(alignment)

        self.assertEqual(expected, actual)

    def test___get_consensus_from_MSA___all_gap_nonmatch(self):
        alignment = make_alignment(["A--A", "A--A"])

        expected = "A**A"
        actual = get_consensus_from_MSA(alignment)

        self.assertEqual(expected, actual)


class TestGetMajorityConsensusFromMSA(TestCase):
    def setUp(self):
        # Create some sample alignments for testing
        self.empty_alignment = MSA([])
        self.all_N_alignment = make_alignment(["NNNN", "NNNN"], ["seq1", "seq2"])
        self.one_N_one_base_alignment = make_alignment(
            ["NNNN", "ACGT"], ["seq1", "seq2"]
        )
        self.with_some_columns_full_of_gaps_and_Ns = make_alignment(
            ["ANG-", "ANG-", "ANN-", "A---", "C---"],
            ["seq1", "seq2", "seq3", "seq4", "seq5"],
        )
        self.sample_alignment = make_alignment(
            ["ACGT", "ACGT", "A--C", "A--G", "C---"],
            ["seq1", "seq2", "seq3", "seq4", "seq5"],
        )
        self.lower_case_sample_alignment = make_alignment(
            ["acgt", "acgt", "a--c", "a--g", "c---"],
            ["seq1", "seq2", "seq3", "seq4", "seq5"],
        )
        self.equally_likely_bases = make_alignment(
            ["AAAAAAAAAAAAAAAAAAAAAAAAA", "CCCCCCCCCCCCCCCCCCCCCCCCC"],
            ["seq1", "seq2"],
        )

    def test___get_majority_consensus_from_MSA___empty_alignment(self):
        consensus = get_majority_consensus_from_MSA(self.empty_alignment)
        self.assertEqual(consensus, "")

    def test___get_majority_consensus_from_MSA___all_Ns(self):
        consensus = get_majority_consensus_from_MSA(self.all_N_alignment)
        for char in consensus:
            self.assertIn(char, "ACGT")

    def test___get_majority_consensus_from_MSA___one_N_one_base(self):
        consensus = get_majority_consensus_from_MSA(self.one_N_one_base_alignment)
        self.assertEqual(consensus, "ACGT")  # The consensus should be the base sequence

    def test___get_majority_consensus_from_MSA___with_some_columns_full_of_gaps_and_Ns(
        self,
    ):
        consensus = get_majority_consensus_from_MSA(
            self.with_some_columns_full_of_gaps_and_Ns
        )
        self.assertEqual(consensus, "AAGC")  # The consensus should be the base sequence

    def test___get_majority_consensus_from_MSA___sample_alignment(self):
        consensus = get_majority_consensus_from_MSA(self.sample_alignment)
        self.assertEqual(
            consensus, "ACGT"
        )  # The consensus should be the most common base at each position

    def test___get_majority_consensus_from_MSA___lower_case_sample_alignment(self):
        consensus = get_majority_consensus_from_MSA(self.lower_case_sample_alignment)
        self.assertEqual(
            consensus, "ACGT"
        )  # The consensus should be the most common base at each position

    def test___get_majority_consensus_from_MSA___equally_likely_bases(self):
        consensus = get_majority_consensus_from_MSA(self.equally_likely_bases)
        self.assertEqual(
            set(consensus), {"A", "C"}
        )  # A and C should be chosen at random
