# Error raising class
class BadSequenceException(Exception):
    pass

# Function for reading a DNA sequence from a FASTA file
def readDNAsequence (file):
    """
    This function takes as its argument the name of a fasta file. 
    When passed the name of a FASTA file, the function reads the file, 
    discards the header and returns the sequence as a string.
    It raises an error if the sequence part of the file contains 
    characters that are not one of the letters A, C, T, G, U. 
    All U nucleotides are replaced by T in the returned string.
    """

    #Importing, calling and reading fasta file
    from IPython.display import FileLink
    FileLink(file)

    file_opened=open(file, "r")

    sequence=file_opened.readlines()

    #Remove header
    sequence.remove(sequence[0])
    sequence

    #Remove \n
    cleaned_sequence = [s.strip() for s in sequence]
    cleaned_sequence

    #String together the cleaned sequence
    stringed_sequence = ''.join(cleaned_sequence)
    stringed_sequence

    #list the sequence
    listed_string = list(stringed_sequence)
    listed_string

    #Account for non-standard nucleotides
    standard_nucleotides = ['A', 'T', 'C', 'G', 'U']
    for i in listed_string:
        if i not in standard_nucleotides:
            raise BadSequenceException("Error: Non-standard characters in sequence or multiple sequences in file.")
    
    #Accounting for Uracils
    if 'U' in listed_string:
        listed_string = [x.replace('U', 'T') for x in listed_string]
        return (''.join(listed_string))
    else:
        return (''.join(listed_string))
    
# Function for computing the complement of a DNA sequence
def complement (sequence):
    """
    This function takes a string containing a DNA sequence as its only 
    parameter and returns the complement of the sequence in a string. 
    The function should raises an error if the argument sequence contains 
    anything else than the four characters A, C, T, G.
    """

    #Turn string sequence into list to be indexed
    listed_sequence = list (sequence)

    #Build nucleotide dictionary
    nucleotide_dictonary = {}
    nucleotide_dictonary ['A'] = 'T'
    nucleotide_dictonary ['T'] = 'A'
    nucleotide_dictonary ['G'] = 'C'
    nucleotide_dictonary ['C'] = 'G'

    #Index nucleotide dictionary for matching bases in the listed sequence
    comp_sequence = [nucleotide_dictonary[x] for x in listed_sequence]

    #Account for non-standard DNA nucleotides
    standard_dna_nucleotides = ['A', 'T', 'C', 'G']
    for i in comp_sequence:
        if i not in standard_dna_nucleotides:
            raise BadSequenceException("Non-standard characters in DNA sequence")
    
    #Print outcome
    return (''.join(comp_sequence))

# Function for extracting primers
def primer (sequence, forward = True, length = 20):
    """
    This function takes three parameters: a DNA sequence (sequence), 
    an integer (length) that is 20 by default and a Boolean value 
    (forward) that is TRUE by default. When forward is set to TRUE 
    (or is not passed), the function returns a Forward primer for the 
    argument sequence. When forward it set to FALSE, it returns a 
    Reverse primer. The length of the primer is specified by the length
    argument. If this is not passed, a primer of length 20 is returned. 
    If the sequence is shorter than the length of nucleotides specified, 
    an error will be raised.
    """

    #Build nucleotide dictionary
    nucleotide_dictionary = {}
    nucleotide_dictionary ['A'] = 'T'
    nucleotide_dictionary ['T'] = 'A'
    nucleotide_dictionary ['G'] = 'C'
    nucleotide_dictionary ['C'] = 'G'

    #if statement to account for forward argument
    if forward:
        primer_sequence = sequence[:length] #Obtaining first 20 characters of string
        
    else:
        raw_reverse_primer_sequence = sequence[-length:] #Obtaining last 20 characters of string
        listed_reverse_primer_sequence = list(raw_reverse_primer_sequence) #listing last 20 elements of string to enable indexing
        listed_reverse_primer_sequence.reverse() #indexed list is reversed to obtain reverse sequence
        reverse_comp_sequence = [nucleotide_dictionary[x] for x in listed_reverse_primer_sequence] #Index dictionary to obtain matching bases
        primer_sequence = ''.join(reverse_comp_sequence) #prints reverse primer sequence

    #Accounting for primer length being shorter than input sequence length
    if len(primer_sequence) < length:
        raise BadSequenceException("Primer length not satisfied")
    
    #Returning output
    return(primer_sequence)

# Function for computing melting temperature
def meltingTemp (primer):
    """
    This function takes a string representing a primer as its argument. 
    The function returns the melting temperature of the primer in degrees Celsius. 
    If the sequence contains characters other than A, C, T, G, the function should raise an error.
    """

    #Account for non-standard DNA nucleotides
    standard_dna_nucleotides = ['A', 'T', 'C', 'G']
    for i in primer:
        if i not in standard_dna_nucleotides:
            raise BadSequenceException("Non-standard characters in DNA sequence")
    
    #Build nucleotide dictionary and obtain total count of each base
    nucleotide_dictionary = {}
    nucleotide_dictionary ['A'] = primer.count("A")
    nucleotide_dictionary ['T'] = primer.count("T")
    nucleotide_dictionary ['G'] = primer.count("G")
    nucleotide_dictionary ['C'] = primer.count("C")
    
    #Calculate melting temperature
    melting_temp = ((4*nucleotide_dictionary['G']) + (4*nucleotide_dictionary['C'])) + ((2*nucleotide_dictionary['A']) + (2*nucleotide_dictionary['T']))

    #Return outcome
    return (melting_temp)

# Function for average melting temperature
def sequencePCRtemp (file):
    """
    This function combines the previous functions.
    It takes a string containing the name of a FASTA file as its argument. 
    It returns the average melting temperature of the two primers of the sequence as a float.
    """

    #Generating Pure DNA Nucleotide Sequence
    pure_dna_sequence = readDNAsequence (file)

    #Primer Extraction
    forward_primer = primer (pure_dna_sequence, forward = True)
    reverse_primer = primer (pure_dna_sequence, forward = False)
    
    #Primer Melting Temperature Calculations
    melting_temp_forward_primer = meltingTemp (forward_primer)
    melting_temp_reverse_primer = meltingTemp (reverse_primer)

    #Average Melting Temperature Calculation
    avg_melting_temp = melting_temp_forward_primer + melting_temp_reverse_primer / 2
    avg_melting_temp = float (avg_melting_temp)
    
    #Return Average Primer Melting Temperature
    return (f"{avg_melting_temp} °C")