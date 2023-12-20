import argparse
import pandas
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from HLA_util import find_amino_acids, find_all_amino_acids_fast

parser = argparse.ArgumentParser(
    description=''' 
    Tool to analyse the amino-acid composition of aligned IPD-IMGT/HLA files.
    '''
)
parser.add_argument(
    '--input', '-i', nargs=1, required=True,
    help='File containing the aligned amino-acid data in IPD-IMGT/HLA format'
)
parser.add_argument(
    '--position', '-p', nargs=1, required=False,
    help='if set the tool will output detailed analysis of specified amino-acid position'
)
parser.add_argument(
    '--count', '-c', action='store_true',
    help='if set the tool will output the total amount of amino-acids found in file '
         'as well as an histogram of different amino-acids found per position'
)


if __name__ == '__main__':
    args = parser.parse_args()
    file = args.input[0]

    if args.position:
        amino_acids = find_amino_acids(file, int(args.position[0]))
        amino_acids_counts = Counter(amino_acids)
        df = pandas.DataFrame.from_dict(amino_acids_counts, orient='index')
        df.plot(kind='bar')
        plt.show()

    if args.count:
        amino_acids_counter = find_all_amino_acids_fast(file)
        total_amino_acids = Counter()
        variations = []
        for c in amino_acids_counter:
            total_amino_acids.update(c)
            common_variations = 0
            for aa, count in c.items():
                if count > 96:
                    common_variations += 1
            variations.append(common_variations)

        # actual position different from position in variations due to negative start position
        variations = np.array(variations)
        print('max variations = ' + str(np.max(variations)))
        print('average variations = ' + str(np.average(variations)))
