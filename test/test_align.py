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
    expected = np.array([[0., -10., -11., -12., -13.],
                        [-10.,  5., -11., -11., -13.],
                        [-11., -10., 4., -1., -6.],
                        [-12., -12., -8., 5., 4.]], dtype=object)
    # assert that the calculated matrix is as expected
    assert (test._align_matrix == expected).all() == True
    

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
    assert test.align(seq3, seq4) == (17.0, 'MAVHQLIRRP', 'M---QLIRHP')





