# HLA-tool
This tool can analyze amino acid frequency in alignment files of a given position of the reference sequence.

## Dependencies
<li>NumPy</li>
<li>Matplotlib</li>

## Usage
For basic usage of finding amino acid frequency of a specific position like 116:
`python HLA.py -i Data/B_prot.txt -p 116`

### Input
This tool uses a alignment file in IPD-IMGT/HLA 3.54.0 format as input.

### Optional Arguments
`--position`,`-p`: if set analyzes the amino acid frequency of specified position of the reference sequence.

`--count`,`-c`: if set the tool will output a histogram of common variations found per position.

`--cluster`,`-cl`: if set the amino-acid histogram will be clustered together based on a shared property. 
                    Currently class and polarity are the options available.

`--protein_threshold`, `-pt`: if set changes the threshold of minimum amount of proteins 
needed before they are considered expressed by the specific amino acid at the specified location.

`--threshold_test`, `-tt`: if set a plot will be generated showing the amount of protein_ids includes 
for each set threshold up to 200,as well as the amount of duplicate ids shared between amino acids.