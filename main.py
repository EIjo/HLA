import re
import pandas
from collections import Counter
import matplotlib.pyplot as plt
import numpy as np


# inconsistent index naming scheme
def find_amino_acids(file_path, position):
    amino_acids = []
    seq_count = 0
    seq_index = 0
    position_index = 0
    line_start = 0
    reference = []
    index_dict = []
    with open('B_prot.txt', 'r') as file:
        for line in file:
            # get line start information as well as the index of the position
            if line.__contains__('Prot'):
                line_start = line.index(line.split()[1][0])
                line_start_index = int(line.split()[1])
                seq_index = 0
                seq_count = seq_index
                next(file)
            elif re.search('[A-Za-z]+\*[0-9]+:[0-9]+.*', line):
                if seq_index == 0:
                    # can be nicer I think
                    seq = line[line_start:]
                    for i in seq:
                        if i == ' ':
                            line_start += 1
                            line_start_index += 1
                        else:
                            index_dict.append(line_start_index)
                            seq = line[line_start:]
                            break
                    ref = seq.strip().replace(' ', '').replace('.', '')
                    reference.append(ref)

                if line_start_index < position < line_start_index + len(ref):
                    if position_index == 0:
                        position_index = line_start

                        i = 0
                        ## with or without +1??? without more variation but counting looks like +1 is valid??
                        while i < position - line_start_index:
                            a = line[position_index:]
                            if not (line[position_index] == ' ' or line[position_index] == '.'):
                                i += 1
                            position_index += 1

                        # position_index = position - line_start_index + 1
                        # print(line[position_index])
                        amino_acids.append(line[position_index])
                        # print(len(line))
                    else:
                        if position_index < len(line):
                            if line[position_index] == '-':
                                amino_acids.append(amino_acids[0])
                            elif line[position_index] == '\n':
                                amino_acids.append('.')
                                # print(len(line))
                            elif line[position_index] == ' ':
                                amino_acids.append('.')
                                # print(len(line))
                            else:
                                amino_acids.append(line[position_index])
                        else:
                            amino_acids.append('.')
                elif position < line_start_index:
                    break
                seq_index += 1
    return amino_acids

# compare variance positions
# variations = []
# for i in range(0, 347):
#     amino_acids = find_amino_acids('B_prot.txt', i)
#     amino_acids_counts = Counter(amino_acids)
#     common_variations = 0
#     for aa, count in amino_acids_counts.items():
#         if count > 96:
#             common_variations += 1
#     variations.append(common_variations)


# variations = np.array(variations)

# histogram amino acids 116
amino_acids = find_amino_acids('B_prot.txt', 116)
amino_acids_counts = Counter(amino_acids)
df = pandas.DataFrame.from_dict(amino_acids_counts, orient='index')
df.plot(kind='bar')
plt.show()
# print(reference)
