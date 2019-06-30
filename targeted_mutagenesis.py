from itertools import product, combinations
from collections import Counter
import pprint

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
valid_aa = 'GAVLIMFWPSTCYNQDEKRH*'



"""
Get a list of all nucleotides for the amino acids at pos1, pos2, and pos3. Convert to sets
Create a dict of sets for position one, position two, and position 3
Items in a set are unique. We can not have a set that contains two items that are equal. This is important when creating.
use is subset
Find the expanded code that intersects with each list, ie all of the letters contained

"""

####################################################################
# get nucleotide lists

test_set = {'A', 'I', 'V'}


#assert statements checking to make sure the amino acid exixts in the translatable
pos1_set = {nuc[0] for nuc, amino in translation_table.items() if amino in test_set}
pos2_set = {nuc[1] for nuc, amino in translation_table.items() if amino in test_set}
pos3_set = {nuc[2] for nuc, amino in translation_table.items() if amino in test_set}

print('Pos 1:')
print(pos1_set)
print('Pos 2:')
print(pos2_set)
print('Pos 3:')
print(pos3_set)

####################################################################

pos_dict = {
    'pos1' : pos1_set,
    'pos2' : pos2_set,
    'pos3' : pos3_set
}

# This collects all the expanded codes containing all the amino acids
pos1_dict = {key:value for (key,value) in expanded_code.items() if pos_dict['pos1'].issubset(set(expanded_code[key]))}
pos2_dict = {key:value for (key,value) in expanded_code.items() if pos_dict['pos2'].issubset(set(expanded_code[key]))}
# pos3_dict = {key:value for (key,value) in expanded_code.items() if pos_dict['pos3'].issubset(set(expanded_code[key]))}

# print(pos1_code_dict)


# this way will keep the key value pair if any nucleotide in that postion is present
# pos1_dict = {key:value for (key,value) in expanded_code.items() if len(pos_dict['pos1'].intersection(set(expanded_code[key]))) > 0}
# pos2_dict = {key:value for (key,value) in expanded_code.items() if len(pos_dict['pos1'].intersection(set(expanded_code[key]))) > 0}


pos3_dict = {key:value for (key,value) in expanded_code.items() if len(pos_dict['pos3'].intersection(set(expanded_code[key]))) > 0}

comb = (set(product(pos1_dict.keys(), pos2_dict.keys(),pos3_dict.keys())))

for combinations in comb:
    print((combinations))

####################################################################


def get_codon_for_amino_acids(amino_acids):
    """
    :param amino_acids: set
        the amino acids we want to code for, i.e. {'A','I','V'}
    :rtype: set, float
        returns two values the set of most efficient codons for the input set list, e.g. {'RYA', 'RYH', 'RYC', 'RYW', 'RYM', 'RYY', 'RYT'} and the achieved efficiency e.g. 0.75
    """

    pass


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
    # assert get_codon_for_amino_acids({'A', 'I', 'V'}) == ({'RYA', 'RYH', 'RYC', 'RYW', 'RYM', 'RYY', 'RYT'}, 0.75)
    # assert get_codon_for_amino_acids({'M', 'F'}) == ({'WTS', 'WTK', "WTB"}, 0.5)

    # # "frozenset" here since this seems to be the only way to get a set of sets - see https://stackoverflow.com/questions/5931291/how-can-i-create-a-set-of-sets-in-python
    # assert truncate_list_of_amino_acids({'A', 'V', 'I'}) == {frozenset({'V', 'A'}), frozenset({'V', 'I'})}

    pass
