import argparse
import pandas
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
import HLA_util

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
parser.add_argument(
    '--cluster', '-cl', nargs=1, required=False,
    help='if set the amino-acid histogram will be clustered together based on a shared property. '
         'Currently class and polarity are the options available'
)
parser.add_argument(
    '--protein_threshold', '-pt', nargs=1, required=False,
    help='if set changes the threshold of minimum amount of proteins needed '
         'before they are considered expressed by the spesific amino acid at the specified location'
)
parser.add_argument(
    '--threshold_test', '-tt', action='store_true',
    help='if set a plot will be generated showing the amount of protein_ids includes for each set threshold up to 200,'
         'as well as the amount of duplicate ids shared between amino acids'
)

if __name__ == '__main__':
    args = parser.parse_args()
    file = args.input[0]

    if args.position:
        amino_acids, amino_acids_ids = HLA_util.find_amino_acids(file, int(args.position[0]))
        amino_acids_counts = Counter(amino_acids)

        # sort by amino acid counts (by highest count first)
        amino_acids_counts = dict(sorted(amino_acids_counts.items(), key=lambda x: -x[1]))
        amino_acids_ids = HLA_util.sort_by_other_dict(amino_acids_ids, amino_acids_counts)

        if args.protein_threshold:
            HLA_util.get_unique_proteins(amino_acids_ids, int(args.protein_threshold[0]))
        else:
            HLA_util.get_unique_proteins(amino_acids_ids)

        if args.threshold_test:
            t = np.arange(0, 200, 1)
            completeness_list, duplicates_list = HLA_util.calc_threshold_proteins(amino_acids_ids)
            fig, ax1 = plt.subplots()

            color = 'tab:red'
            ax1.set_xlabel('threshold value')
            ax1.set_ylabel('%completeness', color=color)
            ax1.plot(t, completeness_list, color=color)
            ax1.tick_params(axis='y', labelcolor=color)

            ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

            color = 'tab:blue'
            ax2.set_ylabel('#duplicates', color=color)  # we already handled the x-label with ax1
            ax2.plot(t, duplicates_list, color=color)
            ax2.tick_params(axis='y', labelcolor=color)

            fig.tight_layout()  # otherwise the right y-label is slightly clipped
            plt.show()

        if args.cluster:
            amino_acids_counts_cluster = HLA_util.cluster_counter(amino_acids_counts, args.cluster[0])

            # sort by amino acid cluster counts (by highest count first)
            amino_acids_counts = dict(sorted(amino_acids_counts.items(), key=lambda x: -x[1]))

            df = pandas.DataFrame.from_dict(amino_acids_counts_cluster, orient='index')
            df.plot(kind='bar')
            plt.show()
        else:
            df = pandas.DataFrame.from_dict(amino_acids_counts, orient='index')
            df.plot(kind='bar')
            plt.show()

    if args.count:
        amino_acids_counter = HLA_util.find_all_amino_acids_fast(file)

        total_amino_acids = Counter()
        variations = []
        for c in amino_acids_counter:
            total_amino_acids.update(c)
            common_variations = 0
            for aa, count in c.items():
                # common variation needs to be higher than 1% of all records
                if count > int(amino_acids_counter[0].total()/100):
                    common_variations += 1
            variations.append(common_variations)

        # actual position different from position in variations due to negative start position
        variations = np.array(variations)
        print('max variations = ' + str(np.max(variations)))
        print('average variations = ' + str(np.average(variations)))



