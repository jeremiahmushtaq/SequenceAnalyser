# Import libraries
import pandas as pd

# Define error raising class
class DimensionalityException(Exception):
    pass

# Function for reading an amino acid sequence
def readAAsequence (file):
    """
    This function takes as its argument the name of a file.  
    When passed the name of a FASTA file containing an aminoacid sequence, 
    the function reads the file, discards the header and returns the sequence as a string.
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
    return (''.join(listed_string))

# Function for amino acid usage statistics
def AAtypes (aa_seq):
    """
    This function takes as its argument a string containing a sequence of aminoacids. 
    The function computes, for each protein, the fraction of aminoacids that are polar, 
    small and hydrophobic and return that as a list.
    """
    
    #Build list of relavent residues
    polar = 'HKRWYDENQSTC'
    small = 'PNDVTSCGA'
    hydrophobic = 'ILVFYWHKMCTA'

    #Segregate relavent residues in amino acid sequence
    polar_in_seq = [x for x in aa_seq if x in polar]
    small_in_seq = [x for x in aa_seq if x in small]
    hydrophobic_in_seq = [x for x in aa_seq if x in hydrophobic]

    #Calculating fractions
    fraction_polar = len(polar_in_seq) / len(aa_seq)
    fraction_small = len(small_in_seq) / len(aa_seq)
    fraction_hydrophobic = len(hydrophobic_in_seq) / len(aa_seq)

    #Assemble in list and limit to two significant figures
    fractions = [round(fraction_polar, 2), round(fraction_small, 2), round(fraction_hydrophobic, 2)]

    #Return fractions
    return (fractions)


# Function for computing Euclidean distance
def distance (a,b):
    """
    This function takes two lists of float type that are of equal length as its arguments. 
    The function returns the Euclidean distance between the two lists, considered as vectors 
    in a space of dimensionality equal to the length of the lists. That is to say, computes 
    the element-by-element difference between the two lists, squares that, sum all the squares 
    and takes a square root of the final sum. The function raises an error if the two lists passed 
    are not of the same length, or if they are both empty.
    """

    #Accounting for unequal input list elements
    if len(a) != len(b):
        raise DimensionalityException("Lists lengths are not equal")

    #Accounting for empty input lists
    if a == [] and b == []:
        raise DimensionalityException("Lists are empty")

    #Calculate difference between each element of lists
    difference = []
    for i, j in zip(a, b):
        x = (i-j)
        difference.append(x)
    
    #Calculate squares
    square = []
    for i in difference:
        x = i**2
        square.append(x)

    #Calculate square root of the sum of squares
    answer = (sum(square))**0.5
    answer = float("%.2f" % answer)
    return(answer)

# Function for computing distance matrix
def DistanceMatrix (filelist):
    """
    This function takes a list of strings called 'filelist' and a separate string called 'outputfile'. 
    Each element of 'filelist' represents the name of a FASTA file containing an aminoacid sequence. 
    For each filename, the function loads the sequence, computes the fraction of aminoacids that are 
    polar, small and/or hydrophobic, and finally output this to the text file specified by 'outputfile' 
    in TSV (tab-separated values) tabular format. The table contains, on each line, the name of the input 
    file and then the percentages of polar, small and hydrophobic residues (in that order).
    """
    
    #Loading and storing amino acid sequences
    aa_seqs = []
    for i in filelist:
        x = readAAsequence(i)
        aa_seqs.append(x)

    #Obtaining fractions from each sequence
    aa_fractions = []
    for i in aa_seqs:
        x = AAtypes(i)
        aa_fractions.append(x)

    #Correspond file names to sorted list
    dict = {}
    for name, sorted in zip(filelist, aa_fractions): 
        dict[name] = sorted

    #Extract and list percentage values from dictionary 
    percentages = list(dict.values())

    #Obtain number of elements in list
    number_of_elements = len(percentages)

    #Calculator distances for each element of the list against each other
    distances = []
    for j in range(0,number_of_elements):
        for i in range(0,number_of_elements):
            x = distance(percentages[j],percentages[i])
            distances.append(x)

    #Sublist distances for each row in the matrix
    distances = [distances[i:i+number_of_elements] for i in range(0, len(distances), number_of_elements)]

    #Extract file names from dictionary
    names = dict.keys()
    names = list(names)

    #Listing file names and corresponding distances
    merged_dict = {name: dist for name, dist in zip(names, distances)}

    #Implement DataFrame
    results_table = pd.DataFrame(merged_dict, index=names).T

    #Convert to HTML format
    html = results_table.to_html()
    
    return(html)