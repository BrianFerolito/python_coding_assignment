from itertools import product, combinations
from collections import Counter
from os.path import join


translation_table = {'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',
                     'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C',
                     'TTA': 'L', 'TCA': 'S', 'TAA': '*', 'TGA': '*',
                     'TTG': 'L', 'TCG': 'S', 'TAG': '*', 'TGG': 'W',
                     'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R',
                     'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',
                     'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',
                     'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',
                     'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S',
                     'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',
                     'ATA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',
                     'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',
                     'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G',
                     'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',
                     'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',
                     'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'
                     }

# nomenclature for degenerate codons
expanded_code = {'A': ['A'], 'C': ['C'], 'G': ['G'], 'T': ['T'],
                 'W': ['A', 'T'], 'S': ['C', 'G'], 'M': ['A', 'C'], 'K': ['G', 'T'], 'R': ['A', 'G'], 'Y': ['C', 'T'],
                 'B': ['C', 'G', 'T'], 'D': ['A', 'G', 'T'], 'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G'],
                 'N': ['A', 'C', 'G', 'T']
                 }

# helpful for validating input
valid_nucleotides = 'ACGTWSMKRYBDHVN'
valid_nucleotides_set = set(list(valid_nucleotides))
valid_aa = 'GAVLIMFWPSTCYNQDEKRH*'
valid_aa_set = set(list(valid_aa))

def get_codon_for_amino_acids(amino_acids):
    """
    :param amino_acids: set
        the amino acids we want to code for, i.e. {'A','I','V'}
    :rtype: set, float
        returns two values the set of most efficient codons for the input set list, e.g. {'RYA', 'RYH', 'RYC', 'RYW', 'RYM', 'RYY', 'RYT'} and the achieved efficiency e.g. 0.75
    """
    
    assert amino_acids.issubset(valid_aa_set), 'Invalid amino acid passed to get_codon_for_amino_acids()'

    # number of amino acids in the set
    amino_length = len(amino_acids)

    # obtain all nucleotides possible at each position
    pos1_set = {nuc[0] for nuc, amino in translation_table.items() if amino in amino_acids}
    pos2_set = {nuc[1] for nuc, amino in translation_table.items() if amino in amino_acids}
    pos3_set = {nuc[2] for nuc, amino in translation_table.items() if amino in amino_acids}

    # narrow done the possible key value pairs from the expanded code dict for each nucleotide position
    pos1_dict = {key:value for (key,value) in expanded_code.items() if set(expanded_code[key]).issubset(pos1_set)}
    pos2_dict = {key:value for (key,value) in expanded_code.items() if set(expanded_code[key]).issubset(pos2_set)}
    pos3_dict = {key:value for (key,value) in expanded_code.items() if set(expanded_code[key]).issubset(pos3_set)}
    
    # get all possible degenerate codons and store the combinations in a list
    comb = (set(product(pos1_dict.keys(), pos2_dict.keys(),pos3_dict.keys())))
    comb_list = [combinations for combinations in comb]

    encoded_dict = {}

    # Loop through every possible degenerate codon
    for degen in comb_list:
        
        # obtain sets of all triplets created by the degenerate codon
        trip_set = (set(product(expanded_code[degen[0]], expanded_code[degen[1]], expanded_code[degen[2]])))
        trip_strings = [''.join(triplets) for triplets in trip_set]
        
        # get set of amino acids from the triplets
        amino_acid_set = {amino for trip, amino in translation_table.items() if trip in trip_strings}

        # if there is a stop codon produced by the degenerate codon skip over it
        if '*' not in amino_acids and '*' in amino_acid_set:
            continue

        # make sure all the amino acids we want to be coded are present in set
        if amino_acids.issubset(amino_acid_set): 
            # add the degenerate codon as the key with the efficiency as the value
            encoded_dict[''.join(degen)] = round(amino_length/len(amino_acid_set), 2)

    # sort the dict by value and return a list of tuples
    ordered_list = sorted(encoded_dict.items(), key=lambda x: x[1], reverse=True)

    # obtain the highest efficiency from the first element of the descending ordered list
    highest_efficiency = ordered_list[0][1]

    # obtain all degenerate codons with the highest efficiency
    efficient_codons = {codon[0] for codon in ordered_list if codon[1] == highest_efficiency}

    return((efficient_codons, highest_efficiency))


def truncate_list_of_amino_acids(amino_acids):
    """
    :param amino_acids: set
        the amino acids we want to code for, i.e. {'A','I','V'}
    :rtype: set
        the set of sets of amino acids that can be coded with 100% efficiency, i.e. {frozenset({'V', 'A'}), frozenset({'V', 'I'})}
    """

    pass


if __name__ == "__main__":
    # using sets instead of lists throughout the code since the order doesn't matter and all items should be unique
    assert get_codon_for_amino_acids({'A', 'I', 'V'}) == ({'RYA', 'RYH', 'RYC', 'RYW', 'RYM', 'RYY', 'RYT'}, 0.75)
    assert get_codon_for_amino_acids({'M', 'F'}) == ({'WTS', 'WTK', "WTB"}, 0.5)

    print(get_codon_for_amino_acids({'A', 'I', 'V'}))
    print(get_codon_for_amino_acids({'M', 'F'}))
    print(get_codon_for_amino_acids({'L', 'R', 'W', 'Y', '*'}))

    # # "frozenset" here since this seems to be the only way to get a set of sets - see https://stackoverflow.com/questions/5931291/how-can-i-create-a-set-of-sets-in-python
    # assert truncate_list_of_amino_acids({'A', 'V', 'I'}) == {frozenset({'V', 'A'}), frozenset({'V', 'I'})}

 
