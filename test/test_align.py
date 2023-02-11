# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")

    # initialize NeedlemanWunsch with appropriate substitution file and gap opening/extending
    test = NeedlemanWunsch(sub_matrix_file = "./substitution_matrices/BLOSUM62.mat", gap_open=-10, gap_extend=-1)
    # align the sequences
    test.align(seq1, seq2)
    # initialize expected matrix
    expected_align = np.array([[0., -10., -11., -12., -13.],
                                [-10.,  5., -11., -11., -13.],
                                [-11., -10., 4., -1., -6.],
                                [-12., -12., -8., 5., 4.]], dtype=object)

    expected_gapA = np.array([[0., -10., -11., -12., -13.],
                             [-np.inf, -21.,  -6.,  -7.,  -8.],
                             [-np.inf, -22., -17.,  -7.,  -8.],
                             [-np.inf, -23., -18., -18.,  -6.]], dtype=object)

    expected_gapB = np.array([[0., -np.inf, -np.inf, -np.inf, -np.inf],
                             [-10., -21., -22., -23., -24.],
                             [-11.,  -6., -17., -18., -19.],
                             [-12.,  -7.,  -7., -12., -17.]], dtype=object)

    # assert that the calculated matrices are as expected
    assert (test._align_matrix == expected_align).all() == True, "Alignment matrix is not as expected."
    assert (test._gapA_matrix == expected_gapA).all() == True, "Gap A matrix is not as expected."
    assert (test._gapB_matrix == expected_gapB).all() == True, "Gap B matrix is not as expected."
    

def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")
    test = NeedlemanWunsch(sub_matrix_file="./substitution_matrices/BLOSUM62.mat", gap_open=-10, gap_extend=-1)
    # assert if the aligned sequences are as expected
    assert test.align(seq3, seq4) == (17.0, 'MAVHQLIRRP', 'M---QLIRHP'), "Alignment between seq3 and seq4 do not yield expected results"





