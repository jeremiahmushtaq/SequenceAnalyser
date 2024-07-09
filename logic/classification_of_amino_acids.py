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

    #Changing fractions to 2 decimal points
    polar_2dp = ("%.2f" % fraction_polar)
    small_2dp = ("%.2f" % fraction_small)
    hydrophobic_2dp = ("%.2f" % fraction_hydrophobic)

    #Assemble in list
    fractions = [polar_2dp, small_2dp, hydrophobic_2dp]

    #Return fractions
    return (fractions)

# Function for processing multiple amino acid sequences
def AAtypetable (filelist, outputfile):
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
        
    #Obtain information for each coloumn
    col_polar = []
    col_small = []
    col_hydro = []

    for i in aa_fractions:
        x = i[0]
        y = i[1]
        z = i[2]
        col_polar.append(x)
        col_small.append(y)
        col_hydro.append(z)
        
    #Writing TSV file
    OUTF = open (outputfile, "w")
    header = "# " + str(outputfile) + "\t" + "Polar" + "\t" + "Small" + "\t" + "Hydro" + "\n"
    OUTF.write(header)
    for a, x, y, z in zip(filelist, col_polar, col_small, col_hydro):
        content = str(a) + "\t" + str(x) + "\t" + str(y) + "\t" + str(z) + "\n"
        OUTF.write(content)
    OUTF.close()
    
    #Return command
    return([])