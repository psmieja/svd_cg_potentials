
_AA_1L_TO_NAME = {
    'A': 'Alanine',
    'R': 'Arginine',
    'N': 'Asparagine',
    'D': 'Aspartic acid',
    'C': 'Cysteine',
    'Q': 'Glutamine',
    'E': 'Glutamic acid',
    'G': 'Glycine',
    'H': 'Histidine',
    'I': 'Isoleucine',
    'L': 'Leucine',
    'K': 'Lysine',
    'M': 'Methionine',
    'F': 'Phenylalanine',
    'P': 'Proline',
    'O': 'Pyrrolysine',
    'S': 'Serine',
    'U': 'Selenocysteine',
    'T': 'Threonine',
    'W': 'Tryptophan',
    'Y': 'Tyrosine',
    'V': 'Valine',
    'B': 'Aspartic acid or Asparagine',
    'Z': 'Glutamic acid or Glutamine',
    'X': 'Any amino acid',
    'J': 'Leucine or Isoleucine',
    '-': 'termination codon'
}

_AA_3L_TO_NAME = {
    'Ala': 'Alanine',
    'Arg': 'Arginine',
    'Asn': 'Asparagine',
    'Asp': 'Aspartic acid',
    'Cys': 'Cysteine',
    'Gln': 'Glutamine',
    'Glu': 'Glutamic acid',
    'Gly': 'Glycine',
    'His': 'Histidine',
    'Ile': 'Isoleucine',
    'Leu': 'Leucine',
    'Lys': 'Lysine',
    'Met': 'Methionine',
    'Phe': 'Phenylalanine',
    'Pro': 'Proline',
    'Pyl': 'Pyrrolysine',
    'Ser': 'Serine',
    'Sec': 'Selenocysteine',
    'Thr': 'Threonine',
    'Trp': 'Tryptophan',
    'Tyr': 'Tyrosine',
    'Val': 'Valine',
    'Asx': 'Aspartic acid or Asparagine',
    'Glx': 'Glutamic acid or Glutamine',
    'Xaa': 'Any amino acid',
    'Xle': 'Leucine or Isoleucine',
    'TERM': 'termination codon'
}

AMINOACIDS = ['G', 'Q', 'C', 'S', 'V', 'L', 'F', 'P', 'M', 'R', 'N', 'Y', 'E', 'A', 'D', 'H', 'T', 'I', 'K', 'W']

def aminoacid_name_from_1l(code):
    return _AA_1L_TO_NAME[code]

def aminoacid_name_from_3l(code):
    return _AA_3L_TO_NAME[code]