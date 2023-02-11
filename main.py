# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    # TODO Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    NW = NeedlemanWunsch(sub_matrix_file = "./substitution_matrices/BLOSUM62.mat", gap_open = -10, gap_extend = -1)
    # align each species to human
    gg_aligned = NW.align(hs_seq, gg_seq)
    mm_aligned = NW.align(hs_seq, mm_seq)
    br_aligned = NW.align(hs_seq, br_seq)
    tt_aligned = NW.align(hs_seq, tt_seq)
    # put all alignment results into a list
    aligned = [(gg_aligned, 'gallus gallus'), (mm_aligned, 'mus musculus'), (br_aligned, 'balaeniceps rex'), (tt_aligned, 'tursiops truncatus')]
    # sort by the score
    aligned_sorted = sorted(aligned, key = lambda x: x[0][0], reverse=True)
    # print species in order of most similar, highest score
    print('species in order of most similar to human BRD2 sequence')
    print([alignment[1] for alignment in aligned_sorted])

    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    # print the comparison and the score
    print('alignment scores corresponding to each species alignment to the human BRD2 sequence')
    print([(alignment[1], alignment[0][0]) for alignment in aligned_sorted])
    

if __name__ == "__main__":
    main()
