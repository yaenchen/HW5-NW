# Importing Dependencies
import numpy as np
from typing import Tuple

# Defining class for Needleman-Wunsch Algorithm for Global pairwise alignment
class NeedlemanWunsch:
    """ Class for NeedlemanWunsch Alignment

    Parameters:
        sub_matrix_file: str
            Path/filename of substitution matrix
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty

    Attributes:
        seqA_align: str
            seqA alignment
        seqB_align: str
            seqB alignment
        alignment_score: float
            Score of alignment from algorithm
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty
    """
    def __init__(self, sub_matrix_file: str, gap_open: float, gap_extend: float):
        # Init alignment and gap matrices
        self._align_matrix = None
        self._gapA_matrix = None
        self._gapB_matrix = None

        # Init matrices for backtrace procedure
        self._back = None
        self._back_A = None
        self._back_B = None

        # Init alignment_score
        self.alignment_score = 0

        # Init empty alignment attributes
        self.seqA_align = ""
        self.seqB_align = ""

        # Init empty sequences
        self._seqA = ""
        self._seqB = ""

        # Setting gap open and gap extension penalties
        self.gap_open = gap_open
        assert gap_open < 0, "Gap opening penalty must be negative."
        self.gap_extend = gap_extend
        assert gap_extend < 0, "Gap extension penalty must be negative."

        # Generating substitution matrix
        self.sub_dict = self._read_sub_matrix(sub_matrix_file) # substitution dictionary

    def _read_sub_matrix(self, sub_matrix_file):
        """
        DO NOT MODIFY THIS METHOD! IT IS ALREADY COMPLETE!

        This function reads in a scoring matrix from any matrix like file.
        Where there is a line of the residues followed by substitution matrix.
        This file also saves the alphabet list attribute.

        Parameters:
            sub_matrix_file: str
                Name (and associated path if not in current working directory)
                of the matrix file that contains the scoring matrix.

        Returns:
            dict_sub: dict
                Substitution matrix dictionary with tuple of the two residues as
                the key and score as value e.g. {('A', 'A'): 4} or {('A', 'D'): -8}
        """
        with open(sub_matrix_file, 'r') as f:
            dict_sub = {}  # Dictionary for storing scores from sub matrix
            residue_list = []  # For storing residue list
            start = False  # trigger for reading in score values
            res_2 = 0  # used for generating substitution matrix
            # reading file line by line
            for line_num, line in enumerate(f):
                # Reading in residue list
                if '#' not in line.strip() and start is False:
                    residue_list = [k for k in line.strip().upper().split(' ') if k != '']
                    start = True
                # Generating substitution scoring dictionary
                elif start is True and res_2 < len(residue_list):
                    line = [k for k in line.strip().split(' ') if k != '']
                    # reading in line by line to create substitution dictionary
                    assert len(residue_list) == len(line), "Score line should be same length as residue list"
                    for res_1 in range(len(line)):
                        dict_sub[(residue_list[res_1], residue_list[res_2])] = float(line[res_1])
                    res_2 += 1
                elif start is True and res_2 == len(residue_list):
                    break
        return dict_sub

    def _check_match(self, value1, value2):
        """
                This function checks if the values are a match or mismatch
                if match return match score
                if mismatch return mismatch score

                Parameters:
                	value1: str
                 		the first string
                 	value2: str
                 		the second string to be checked

                Returns:
                 	respective score
        """

        if (value1, value2) in self.sub_dict:
            # get the score for the two values, return error if not a match/mismatch
            return self.sub_dict.get((value1, value2))
        else:
            return 'GAP'

    def align(self, seqA: str, seqB: str) -> Tuple[float, str, str]:
        """
        TODO
        
        This function performs global sequence alignment of two strings
        using the Needleman-Wunsch Algorithm
        
        Parameters:
        	seqA: str
         		the first string to be aligned
         	seqB: str
         		the second string to be aligned with seqA
         
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        # Resetting alignment in case method is called more than once
        self.seqA_align = ""
        self.seqB_align = ""

        # Resetting alignment score in case method is called more than once
        self.alignment_score = 0

        # Initializing sequences for use in backtrace method
        self._seqA = seqA
        self._seqB = seqB

        # create variables for size of matrix based on length of the sequences
        self._size_A = len(seqA)
        self._size_B = len(seqB)

        # create matrices for alignment scores, gaps, and backtracing
        # create matrix for alignment scores
        # need to add one because there will be an extra row/column for the gaps
        self._align_matrix = np.ones([(self._size_B + 1), (self._size_A + 1)]) * -np.inf
        # initialize gap matrices
        self._gapA_matrix = np.ones([(self._size_B + 1), (self._size_A + 1)]) * -np.inf
        self._gapB_matrix = np.ones([(self._size_B + 1), (self._size_A + 1)]) * -np.inf

        # initialize matrices for backtrace procedure
        self._back = np.empty([(self._size_B + 1), (self._size_A + 1)])
        self._back_A = np.empty([(self._size_B + 1), (self._size_A + 1)])
        self._back_B = np.empty([(self._size_B + 1), (self._size_A + 1)])

        # in matrix for alignment scores
        # initialize the first row to the gap penalty
        self._align_matrix[0] = np.insert([self.gap_open + self.gap_extend * i for i in range(0, self._size_A)], [0][0], 0)
        # initialize the first column to the gap penalty
        self._align_matrix[:, 0] = np.insert([self.gap_open + self.gap_extend * i for i in range(0, self._size_B)], [0][0], 0)

        # initialize gapA and gapB matrices
        # gapA will hold scores based on the gap opening/extension of sequence A
        self._gapA_matrix[0] = np.insert([self.gap_open + self.gap_extend * i for i in range(0, self._size_A)], [0][0], 0)
        # gapB will hold scores based on the gap opening/extension of sequence B
        self._gapB_matrix[:, 0] = np.insert([self.gap_open + self.gap_extend * i for i in range(0, self._size_B)], [0][0], 0)

        # add all the other values in the matrix
        # take the maximum of the top-left, left, and top neighbors
        for row in range(1, self._size_B + 1):
            for column in range(1, self._size_A + 1):
                # get the maximum values for the alignment and gap matrices
                # get the options first, starting from top left
                align_options = [self._align_matrix[row - 1,column - 1],
                                 self._gapA_matrix[row - 1, column - 1],
                                 self._gapB_matrix[row - 1, column - 1]]
                # get the maximum value
                max_align = max(align_options)
                # repeat for gap matrices
                a_options = [self.gap_open + self.gap_extend + self._align_matrix[row, column - 1],
                             self.gap_extend + self._gapA_matrix[row, column - 1],
                             self.gap_open + self.gap_extend + self._gapB_matrix[row, column - 1]]
                max_a = max(a_options)
                b_options = [self.gap_open + self.gap_extend + self._align_matrix[row - 1, column],
                             self.gap_open + self.gap_extend + self._gapA_matrix[row - 1, column],
                             self.gap_extend + self._gapB_matrix[row - 1, column]]
                max_b = max(b_options)

                # calculate values for each matrix using the maximum values
                # for alignment matrix, also add the mismatch/match value from the substitution BLOSUM/PAM matrix
                self._align_matrix[row, column] = self.sub_dict.get((self._seqB[row - 1], self._seqA[column - 1])) + max_align
                # for gap matrices, first case is always a match, second case is a gap in seqA, last case is a gap in seqB
                self._gapA_matrix[row, column] = max_a
                self._gapB_matrix[row, column] = max_b

                # update the backtracking matrices
                self._back[row, column] = align_options.index(max_align)
                self._back_A[row, column] = a_options.index(max_a)
                self._back_B[row, column] = b_options.index(max_b)

                # # check if match, mismatch, or gap
                # new_score = self._check_match(self._seqB[row - 1], self._seqA[column - 1])
                # # if match or mismatch, add the score from the BLOSUM matrix
                # # if gap, add gap score
                # if new_score == 'GAP':
                #     new_score = self.gap_open
                # # calculate top and left scores, which both represent gaps
                # top = scores[row - 1][column] + self.gap_open
                # left = scores[row][column - 1] + self.gap_open
                # # calculate diagonal top left score based on if it is a match, mismatch, or gap
                # top_left = scores[row - 1][column - 1] + new_score
                #
                # # add the maximum score to the scores matrix
                # scores[row][column] = max(top, left, top_left)
                # # add the specified direction in the direction matrix, which will be used to backtrack
                # if max(top, left, top_left) == top:
                #     direction_matrix[row][column] = 1
                # elif max(top, left, top_left) == left:
                #     direction_matrix[row][column] = 2
                # elif max(top, left, top_left) == top_left:
                #     direction_matrix[row][column] = 3

        #self.score_matrix = scores
        #self.direction_matrix = direction_matrix

        print(self._align_matrix)
        return self._backtrace()

    def _backtrace(self) -> Tuple[float, str, str]:
        """
        TODO
        
        This function traces back through the back matrix created with the
        align function in order to return the final alignment score and strings.
        
        Parameters:
        	None
        
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """

        # initialize starting position
        j = self._size_A # size of sequence, 'ABCD' would be 4
        i = self._size_B
        # get the final alignment score, the bottom right value of the final matrix
        self.alignment_score = max(self._align_matrix[i, j], self._gapA_matrix[i, j], self._gapB_matrix[i, j])

        # j is num col/first sequence, and i is num row/second sequence
        while i > 0 and j > 0:
            # get maximum at each position until you hit top left position where i and j are both 0
            # choose between gaps or match/mismatch
            final_options = [self._align_matrix[i, j],
                             self._gapA_matrix[i, j],
                             self._gapB_matrix[i, j]]
            final_max = max(final_options)
            # keep track of the option so we can backtrace with index
            max_position = final_options.index(final_max)

            # backtrace using the position of the maximum from final_options
            # top left case
            if max_position == 0:
                # go to the top left value
                j -= 1
                i -= 1
                # top left is a match, add same character to both sequences
                self.seqA_align += self._seqA[j]
                self.seqB_align += self._seqB[i]
            # left case
            elif max_position == 1:
                # go to the left value
                j -= 1
                # top is a gap, add a gap to the second sequence
                self.seqB_align += '-'
                self.seqA_align += self._seqA[j]
            # top case
            elif max_position == 2:
                # go to the top value
                i -= 1
                # top is a gap, add a gap to the second sequence and keep the first sequence value
                self.seqB_align += '-'
                self.seqA_align += self._seqA[i + 1]

        # when we reach the most top left cell
        # reverse the final sequences since we worked backwards
        self.seqA_align = self.seqA_align[::-1]
        self.seqB_align = self.seqB_align[::-1]

        # # calculate final alignment score
        # alignment_index = 0
        # # while the index is less than the length of the aligned sequence
        # while alignment_index < len(self.seqB_align):
        #     # if the value at that index is a gap
        #     if self._check_match(self.seqA_align[alignment_index], self.seqB_align[alignment_index]) == 'GAP':
        #         # if it is a single gap, add opening gap penalty score
        #         self.alignment_score += self.gap_open
        #         alignment_index += 1
        #         # if the sequence ends on a gap, stop
        #         if alignment_index >= len(self.seqB_align):
        #             break
        #         # while the gap continues
        #         while self.seqB_align[alignment_index] == '-':
        #             # add the extension penalty score
        #             self.alignment_score += self.gap_extend
        #             alignment_index += 1
        #     # if there was no gap
        #     else:
        #         # add the match/mismatch score from the substitution matrix
        #         self.alignment_score += self._check_match(self.seqA_align[alignment_index], self.seqB_align[alignment_index])
        #         alignment_index += 1
        #
        return (self.alignment_score, self.seqA_align, self.seqB_align)


def read_fasta(fasta_file: str) -> Tuple[str, str]:
    """
    DO NOT MODIFY THIS FUNCTION! IT IS ALREADY COMPLETE!

    This function reads in a FASTA file and returns the associated
    string of characters (residues or nucleotides) and the header.
    This function assumes a single protein or nucleotide sequence
    per fasta file and will only read in the first sequence in the
    file if multiple are provided.

    Parameters:
        fasta_file: str
            name (and associated path if not in current working directory)
            of the Fasta file.

    Returns:
        seq: str
            String of characters from FASTA file
        header: str
            Fasta header
    """
    assert fasta_file.endswith(".fa"), "Fasta file must be a fasta file with the suffix .fa"
    with open(fasta_file) as f:
        seq = ""  # initializing sequence
        first_header = True
        for line in f:
            is_header = line.strip().startswith(">")
            # Reading in the first header
            if is_header and first_header:
                header = line.strip()  # reading in fasta header
                first_header = False
            # Reading in the sequence line by line
            elif not is_header:
                seq += line.strip().upper()  # generating full sequence
            # Breaking if more than one header is provided in the fasta file
            elif is_header and not first_header:
                break
    return seq, header
