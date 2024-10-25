"""
bio_functions.py

This module provides core bioinformatics functions for analyzing DNA sequences:
- reverse_complement: Computes the reverse complement of a DNA sequence.
- translate: Translates a DNA sequence into a protein sequence based on the standard genetic code.

Authors:
    - Chris Anderson (jcanderson@berkeley.edu, jcaucb)
    - ChatGPT (OpenAI's GPT-4o)

Implementation Notes:
    This script is designed to provide foundational tools for sequence manipulation in bioinformatics.
    Each function includes validation, error handling, and adheres to biological conventions for sequence analysis.

    - reverse_complement: Computes the reverse complement by reversing the DNA sequence and substituting each nucleotide with its complementary pair.
    - translate: Converts DNA sequences into amino acid sequences using a standard codon table, and handles stop codons by translating them to underscores.

Usage:
    The functions can be imported for use in a larger bioinformatics pipeline or executed directly in this script for demonstration purposes.
"""

def reverse_complement(sequence):
    """
    Calculates the reverse complement of a DNA sequence.

    This function first validates that all characters in the sequence are valid nucleotides (A, T, C, G).
    It then reverses the sequence and replaces each nucleotide with its complementary base:
    A <-> T, C <-> G. This approach ensures the reverse strand matches biological conventions
    for base pairing in DNA.

    Args:
        sequence (str): A string representing the DNA sequence.

    Returns:
        str: The reverse complement of the DNA sequence.

    Raises:
        ValueError: If the DNA sequence contains invalid characters (anything other than A, T, C, G).

    Example:
        reverse_complement("ATGC")  # Returns "GCAT"

    Implementation Notes:
        - This function uses a dictionary lookup to perform complementation, which ensures O(n) time complexity,
          where n is the length of the sequence.
        - Error handling is incorporated to raise a ValueError if any characters are found outside the set {A, T, C, G}.
    """
    # Valid nucleotides for DNA; only A, T, C, G are allowed
    valid_nucleotides = {'A', 'T', 'C', 'G'}
    if any(char not in valid_nucleotides for char in sequence):
        raise ValueError("DNA sequence contains invalid characters. Allowed characters: A, T, C, G.")

    # Complement mapping dictionary for DNA bases
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    # Reversing the sequence and substituting with complements
    return ''.join(complement[base] for base in reversed(sequence))


def translate(sequence):
    """
    Translates a DNA sequence into a protein sequence based on the standard genetic code.

    The function reads the DNA sequence in triplets (codons) and maps each codon to its corresponding amino acid.
    Sequences not divisible by three raise an error, as partial codons are biologically invalid.
    Stop codons (TAA, TAG, TGA) are represented by underscores (_) to denote termination points.

    Args:
        sequence (str): A string representing the DNA sequence.

    Returns:
        str: The corresponding protein sequence.

    Raises:
        ValueError: If the DNA sequence contains invalid characters or if the length is not a multiple of three.

    Example:
        translate("ATGGCC")  # Returns "MA"

    Implementation Notes:
        - This function uses a dictionary (`codon_table`) for efficient O(1) lookup of each codon.
        - Stop codons are represented as underscores to signal the end of a protein sequence.
        - Error handling ensures only valid DNA nucleotides and sequence lengths (multiples of three) are processed.
    """
    # Valid nucleotides; only A, T, C, G are allowed
    valid_nucleotides = {'A', 'T', 'C', 'G'}
    if any(char not in valid_nucleotides for char in sequence):
        raise ValueError("DNA sequence contains invalid characters. Allowed characters: A, T, C, G.")
    if len(sequence) % 3 != 0:
        raise ValueError("Length of DNA sequence is not a multiple of three, which is required for translation.")

    # Standard genetic code codon table
    codon_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    
    # Translate the DNA sequence into a protein sequence by iterating over each codon
    protein = ""
    for i in range(0, len(sequence), 3):
        codon = sequence[i:i+3]
        protein += codon_table.get(codon, '_')  # Stop codons or unknown codons map to '_'
    return protein


if __name__ == "__main__":
    # Example DNA sequence for demonstration
    dna_example = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"

    try:
        # Perform reverse complement
        rc_result = reverse_complement(dna_example)
        print(f"Reverse Complement of '{dna_example}': {rc_result}")

        # Perform translation
        translation_result = translate(dna_example)
        print(f"Translation of '{dna_example}': {translation_result}")
    except Exception as e:
        print(f"Error: {str(e)}")
