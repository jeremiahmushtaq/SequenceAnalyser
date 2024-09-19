# Import function to read DNA sequence
from functions.primers_and_melting_temperature import readDNAsequence
import textwrap


# Error raising class
class OpenReadingFrameException(Exception):
    pass

# Function for reading frames
def translate(dna_sequence):
    """
    This function takes a string containing a DNA sequence as its input 
    and outputs a Python dictionary containing the translation of the sequence 
    in all possible reading frames. The keys of the dictionary are f1, f2 and f3 
    for the three forward frames and r1, r2 and r3 for the reverse reading frames.
    The value of each key is the translation of the sequence in the corresponding frame.
    An asterisk(*) represents stop codons.
    """

    ##############################################
    #Generating mRNA sequence
    ##############################################

    #Replacing Ts with Us
    listed_string = [x.replace('T', 'U') for x in dna_sequence]
    mrna_sequence = ''.join(listed_string)

    ###############################################
    # Generating reading frames
    ###############################################

    #Listing mRNA sequence
    listed_mrna= list(mrna_sequence)

    #Indexing the lists for the forward reading frames
    listed_forward_reading_frame_2 = listed_mrna [1:]
    listed_forward_reading_frame_3 = listed_mrna [2:]

    #Stringing forward reading frames
    f1 = mrna_sequence
    f2 = ''.join(listed_forward_reading_frame_2)
    f3 = ''.join(listed_forward_reading_frame_3)

    #Reversing mRNA sequence
    listed_mrna.reverse() 

    #Indexing the lists for the reverse reading frames
    listed_reverse_reading_frame_1 = listed_mrna [:]
    listed_reverse_reading_frame_2 = listed_mrna [1:]
    listed_reverse_reading_frame_3 = listed_mrna [2:]

    #Stringing reverse reading frames
    r1 = ''.join(listed_reverse_reading_frame_1)
    r2 = ''.join(listed_reverse_reading_frame_2)
    r3 = ''.join(listed_reverse_reading_frame_3)

    ###############################################
    # Seperating codons in reading frame
    ###############################################

    #Accumulators of seperated codons
    f1_listed_codons = []
    f2_listed_codons = []
    f3_listed_codons = []
    r1_listed_codons = []
    r2_listed_codons = []
    r3_listed_codons = []

    #Iterative loops for seperating codons
    for i in range(0, len(f1), 3):
        codons = f1[i:i+3]
        f1_listed_codons.append(codons)

    for i in range(0, len(f2), 3):
        codons = f2[i:i+3]
        f2_listed_codons.append(codons)

    for i in range(0, len(f3), 3):
        codons = f3[i:i+3]
        f3_listed_codons.append(codons)

    for i in range(0, len(r1), 3):
        codons = r1[i:i+3]
        r1_listed_codons.append(codons)

    for i in range(0, len(r2), 3):
        codons = r2[i:i+3]
        r2_listed_codons.append(codons)

    for i in range(0, len(r3), 3):
        codons = r3[i:i+3]
        r3_listed_codons.append(codons)

    #Cleaning non-triplet codons
    for i in f1_listed_codons:
        if len(i) != 3:
            f1_listed_codons.remove(i)

    for i in f2_listed_codons:
        if len(i) != 3:
            f2_listed_codons.remove(i)

    for i in f3_listed_codons:
        if len(i) != 3:
            f3_listed_codons.remove(i)

    for i in r1_listed_codons:
        if len(i) != 3:
            r1_listed_codons.remove(i)

    for i in r2_listed_codons:
        if len(i) != 3:
            r2_listed_codons.remove(i)

    for i in r3_listed_codons:
        if len(i) != 3:
            r3_listed_codons.remove(i)

    ###############################################
    # Building codon-amino acid dictionary
    ###############################################

    codon = {
    'UCA': 'S',    # Serine
    'UCC': 'S',    # Serine
    'UCG': 'S',    # Serine
    'UCU': 'S',    # Serine
    'UUC': 'F',    # Phenylalanine
    'UUU': 'F',    # Phenylalanine
    'UUA': 'L',    # Leucine
    'UUG': 'L',    # Leucine
    'UAC': 'Y',    # Tyrosine
    'UAU': 'Y',    # Tyrosine
    'UAA': '*',    # Stop
    'UAG': '*',    # Stop
    'UGC': 'C',    # Cysteine
    'UGU': 'C',    # Cysteine
    'UGA': '*',    # Stop
    'UGG': 'W',    # Tryptophan
    'CUA': 'L',    # Leucine
    'CUC': 'L',    # Leucine
    'CUG': 'L',    # Leucine
    'CUU': 'L',    # Leucine
    'CCA': 'P',    # Proline
    'CCC': 'P',    # Proline
    'CCG': 'P',    # Proline
    'CCU': 'P',    # Proline
    'CAC': 'H',    # Histidine
    'CAU': 'H',    # Histidine
    'CAA': 'Q',    # Glutamine
    'CAG': 'Q',    # Glutamine
    'CGA': 'R',    # Arginine
    'CGC': 'R',    # Arginine
    'CGG': 'R',    # Arginine
    'CGU': 'R',    # Arginine
    'AUA': 'I',    # Isoleucine
    'AUC': 'I',    # Isoleucine
    'AUU': 'I',    # Isoleucine
    'AUG': 'M',    # Methionine
    'ACA': 'T',    # Threonine
    'ACC': 'T',    # Threonine
    'ACG': 'T',    # Threonine
    'ACU': 'T',    # Threonine
    'AAC': 'N',    # Asparagine
    'AAU': 'N',    # Asparagine
    'AAA': 'K',    # Lysine
    'AAG': 'K',    # Lysine
    'AGC': 'S',    # Serine
    'AGU': 'S',    # Serine
    'AGA': 'R',    # Arginine
    'AGG': 'R',    # Arginine
    'GUA': 'V',    # Valine
    'GUC': 'V',    # Valine
    'GUG': 'V',    # Valine
    'GUU': 'V',    # Valine
    'GCA': 'A',    # Alanine
    'GCC': 'A',    # Alanine
    'GCG': 'A',    # Alanine
    'GCU': 'A',    # Alanine
    'GAC': 'D',    # Aspartic acid
    'GAU': 'D',    # Aspartic acid
    'GAA': 'E',    # Glutamic acid
    'GAG': 'E',    # Glutamic acid
    'GGA': 'G',    # Glycine
    'GGC': 'G',    # Glycine
    'GGG': 'G',    # Glycine
    'GGU': 'G'     # Glycine
    }

    ###############################################
    # Generate corresponding amino acids
    ###############################################

    #Accumulators of translated sequence
    translated_sequence_f1 = []
    translated_sequence_f2 = []
    translated_sequence_f3 = []
    translated_sequence_r1 = []
    translated_sequence_r2 = []
    translated_sequence_r3 = []

    #Iterative loops for generating corresponding amino acid sequences
    for i in f1_listed_codons:
        corresponding_translated_sequence = codon.get(i)
        translated_sequence_f1.append(corresponding_translated_sequence)

    for i in f2_listed_codons:
        corresponding_translated_sequence = codon.get(i)
        translated_sequence_f2.append(corresponding_translated_sequence)    

    for i in f3_listed_codons:
        corresponding_translated_sequence = codon.get(i)
        translated_sequence_f3.append(corresponding_translated_sequence)

    for i in r1_listed_codons:
        corresponding_translated_sequence = codon.get(i)
        translated_sequence_r1.append(corresponding_translated_sequence)    

    for i in r2_listed_codons:
        corresponding_translated_sequence = codon.get(i)
        translated_sequence_r2.append(corresponding_translated_sequence)  

    for i in r3_listed_codons:
        corresponding_translated_sequence = codon.get(i)
        translated_sequence_r3.append(corresponding_translated_sequence)  

    ###############################################
    # Append translated sequences
    ###############################################

    #Accumulators of translated sequences
    stringed_translated_sequence_f1 = []
    stringed_translated_sequence_f2 = []
    stringed_translated_sequence_f3 = []
    stringed_translated_sequence_r1 = []
    stringed_translated_sequence_r2 = []
    stringed_translated_sequence_r3 = []

    #Iterative loops for translating sequences
    for i in translated_sequence_f1:
        stringed_translated_sequence_f1.append(i)

    for i in translated_sequence_f2:
        stringed_translated_sequence_f2.append(i)

    for i in translated_sequence_f3:
        stringed_translated_sequence_f3.append(i)

    for i in translated_sequence_r1:
        stringed_translated_sequence_r1.append(i)

    for i in translated_sequence_r2:
        stringed_translated_sequence_r2.append(i)

    for i in translated_sequence_r3:
        stringed_translated_sequence_r3.append(i)

    ###############################################
    # Building final reading frame dictionary
    ###############################################

    #Build dictionary
    amino_acid_dict = {}
    amino_acid_dict ['f1'] = ''.join(stringed_translated_sequence_f1)
    amino_acid_dict ['f2'] = ''.join(stringed_translated_sequence_f2)
    amino_acid_dict ['f3'] = ''.join(stringed_translated_sequence_f3)
    amino_acid_dict ['r1'] = ''.join(stringed_translated_sequence_r1)
    amino_acid_dict ['r2'] = ''.join(stringed_translated_sequence_r2)
    amino_acid_dict ['r3'] = ''.join(stringed_translated_sequence_r3)

    #Return final dictionary
    return (amino_acid_dict)

# Function for locating open reading frame
def openReadingFrame(sequence):
    """
    This function takes a string containing an aminoacid sequence as its argument 
    and returns a string containing the aminoacids between the first Methionine (included) 
    and the first STOP codon that follows it. The stop codon is represented by 
    an asterisk(*), the same as the translate function above. If either the Methionine or 
    the STOP codon are missing, the function returns an empty string.
    """

    if 'M' in sequence and '*' in sequence:
        met_position = sequence.index('M')
        stop_position = sequence.find('*', met_position)
        if stop_position != -1:
            orf = sequence[met_position: stop_position + 1]
            return(orf)
        else:
            return("")
    else:
        return("")

# Function for translating a DNA sequence
def candidateProtein (dna_sequence):
    """
    This function takes a string containing a DNA sequence as its input and outputs 
    the string of aminoacids corresponding to the longest ORF. This function is dependant
    on the translate and openReadingFrame functions above.
    """

    ####################################################################
    #Translating DNA seqeuence for all possible reading frames
    ####################################################################

    translated_sequence = translate (dna_sequence)
    
    ####################################################################
    #Storing amino acid sequences from relavent dictionary into a list
    ####################################################################

    amino_acid_seqs = list(translated_sequence.values())
    
    ####################################################################
    #Extracting and cleaning ORFs
    ####################################################################

    #Accumulator 
    orfs = [] 

    #Iterative loop for extracting ORFs
    for i in amino_acid_seqs: 
        x = openReadingFrame(i)
        orfs.append(x)

    #Removes empty strings
    orfs = [x for x in orfs if x != '']
    
    ####################################################################
    #Error handing if there is no open reading frame
    ####################################################################
    if len(orfs) == 0:
        raise OpenReadingFrameException("No open reading frame found in sequence.")

    ####################################################################
    #Obtaining open reading frame lengths & calculating max length
    ####################################################################

    #Accumulator
    orf_lengths = [] 

    #Iterative loop for counting ORF lengths
    for i in orfs:
        x = len (i)
        orf_lengths.append(x)

    #Calculating longest ORF length
    max_length = max(orf_lengths)

    ####################################################################
    #Corresponding ORFs to their lengths
    ####################################################################

    #Empty dictionary
    orf_length_dict = {}  

    #Combining keys and values and add them to the dictionary
    for orf_lengths, orfs in zip(orf_lengths, orfs): 
        orf_length_dict[orf_lengths] = orfs
        
    #Indexing dictionary for longest ORF
    longest_ORF = orf_length_dict[max_length]

    # Limit line width
    longest_ORF = textwrap.fill(longest_ORF, width=100)

    #Print longest ORF
    return (longest_ORF)

# Function to write a fasta file
def writeFASTA (sequence, description, filename):
    """
    This function takes three string arguments called, in the order, sequence, description and filename. 
    The sequence argument contains an aminoacid sequence. The description argument contains a description 
    (eg name of protein, organism, etc). The filename argument contains a file name. The function creates 
    a file with the name requested, writes to it the description provided as a FASTA header (i.e. starting 
    with the character '>') and writes the sequence to the file. Long sequences are formatted over several lines.
    """

    #Candidate protein
    filename=filename
    description=description
    sequence=sequence

    #Cleaning amino acid sequence
    sequence=sequence.replace(' ','')
    sequence=sequence.replace('\t', '')
    sequence=sequence.replace('\n', '')

    #Opening the file in writing mode
    OUTF=open(filename, "w")
    
    #Generating header
    header="> " + filename + " | " + description+"\n"
    OUTF.write(header)
    OUTF.write(sequence)
    OUTF.close() #Closing file

# Function for finding the longest ORF of a DNA sequence
def maximalORFa (inputfile):
    """
    This function takes as its argument string 'inputfile' containing the name of an input file.
    The function reads a DNA sequence from the input file and returns the candidate protein 
    corresponding to the longest ORF.
    """

    #Generate pure DNA nucleotide sequence
    pure_dna_sequence = readDNAsequence(inputfile)

    #Generate longest ORF
    longest_ORF = candidateProtein(pure_dna_sequence)

    return (longest_ORF)

# Function for writing a fasta file with the longest ORF of a DNA sequence
def maximalORFb (inputfile, outputfile, proteinname):
    """
    This function takes as its argument string 'inputfile' containing the name of an input file, 
    string 'outputfile' with the name of an output file and string 'proteinname' with a description 
    of a candidate protein. The function reads a DNA sequence from the input file and writes the 
    candidate protein corresponding to the longest ORF to the output file, in FASTA format. 
    The string supplied in 'proteinname' provides the header of the FASTA file.
    """

    #Generate pure DNA nucleotide sequence
    pure_dna_sequence = readDNAsequence(inputfile)

    #Generate longest ORF
    longest_ORF = candidateProtein(pure_dna_sequence)

    #Writing output fasta file
    output_fasta_file = writeFASTA (longest_ORF, proteinname, outputfile)