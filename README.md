# SequenceAnalyser

SequenceAnalyser is a bioinformatics tool designed to help researchers analyse DNA and amino acid sequences. This web application provides three main services:

1. **Primer Melting Temperature Calculation**: Takes a DNA sequence in FASTA format and computes the forward and reverse primers. It then uses the Wallace-Ikatura equation to calculate the melting temperatures of both primers and provides an average melting temperature.
2. **Longest Open Reading Frame Determination**: Accepts a DNA sequence in FASTA format, translates it into amino acids, and identifies the longest open reading frame (ORF), starting with a methionine and ending with a stop codon.
3. **Amino Acid Clustering Matrix**: Accepts multiple amino acid sequences in FASTA format and computes a pairwise differentiation matrix to compare protein sequences based on their amino acid composition and properties (polar, small, hydrophobic). This feature is useful for protein similarity analysis.

## How to Use SequenceAnalyser

Launch the application:  
[SequenceAnalyser](https://sequenceanalyser-449372077532.europe-west2.run.app/)

### Primer Melting Temperature
1. Upload a FASTA file containing a DNA sequence using "Choose file".
2. The file must follow these rules:
    - Accepted formats: `.txt`, `.fasta`, `.fas`, `.fa`, `.fna`, `.ffn`, `.faa`, `.mpfa`, `.frn`.
    - One sequence per file with a standard FASTA header (starting with `>`).
    - The DNA sequence should consist of standard bases (A, C, T, G) without gaps.
3. Click "Upload and Calculate".
4. The output will be the average melting temperature of the forward and reverse primers in degrees Celsius.

### Longest Open Reading Frame
1. Upload a FASTA file containing a DNA sequence using "Choose file".
2. The file must follow these rules:
    - Accepted formats: `.txt`, `.fasta`, `.fas`, `.fa`, `.fna`, `.ffn`, `.faa`, `.mpfa`, `.frn`.
    - One sequence per file with a standard FASTA header (starting with `>`).
    - The DNA sequence should consist of standard bases (A, C, T, G) without gaps.
3. Click "Upload and Analyse".
4. The output will be the longest DNA sequence starting with methionine and ending with an asterisk representing the stop codon.

### Protein Clustering
1. Upload multiple FASTA files containing amino acid sequences using "Choose file".
2. The files must follow these rules:
    - Accepted formats: `.txt`, `.fasta`, `.fas`, `.fa`, `.fna`, `.ffn`, `.faa`, `.mpfa`, `.frn`.
    - One sequence per file with a standard FASTA header (starting with `>`).
    - The amino acid sequences should consist of standard one-letter amino acid codes (e.g., AEISQYW) without gaps.
3. Click "Upload and Cluster".
4. The output will be a pairwise differentiation matrix comparing fractions of polar, small, and hydrophobic amino acids for each protein.

## Limitations
Please note the following limitations when using SequenceAnalyser:

1. **Longest ORF Selection**: The application identifies the first open reading frame (ORF) in each reading frame and compares it with others. Multiple ORFs within a reading frame are not considered.
2. **ORF Start Codon**: SequenceAnalyser assumes that the longest ORF begins with methionine (start codon in eukaryotes and archaea) or N-formylmethionine in bacteria, mitochondria, and plastids.
3. **Primer Melting Temperature**: The melting temperature is an average of the forward and reverse primers, so it may not precisely reflect the melting temperature of either primer individually.
4. **File Format**: The application only accepts one sequence per FASTA file. Multiple sequences in a single file will trigger an error.

## Future Enhancements
We are continuously working to improve SequenceAnalyser. In future updates, we aim to introduce the following features:

1. **Handling Protein Tertiary Structure Data**: A future feature could involve handling protein tertiary structure data (if available), allowing clustering based on structural data rather than just sequence data.
2. **Parallel ORF Detection**: The current implementation only considers the first open reading frame in each reading frame. Future updates will include functionality to detect multiple ORFs within a single reading frame.
3. **Optimization for Large Datasets**: As the application scales to handle larger datasets, we will introduce performance optimizations to ensure efficient processing of sequence data.
4. **Improving Primer Melting Temperature Calculation**: Introduce more advanced primer melting temperature calculation methods that account for sequence context, mismatches, and secondary structure formation. This would allow for more precise estimates of the forward and reverse primer melting temperatures.
5. **Handling Multiple ORFs in a Frame**: Extend the current functionality to detect and report multiple ORFs within each reading frame. This would enable users to explore all potential open reading frames in the sequence, not just the first one identified.
6. **Support for Multiple Sequences in One FASTA File**: Add functionality to handle multiple sequences within a single FASTA file, allowing for batch processing and analysis of multiple DNA or protein sequences in one go.
7. **Flexible Start Codon Detection**: Include an option to detect alternative start codons besides methionine (or N-formylmethionine) for users working with diverse organisms. This would make the tool applicable to a broader range of biological systems.

## Contact & Support
For any issues, suggestions, or feedback, please contact or raise an issue on [Issues](https://github.com/jeremiahmushtaq/SequenceAnalyser/issues).
