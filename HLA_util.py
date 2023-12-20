import re
from collections import Counter


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
                # add reference at start of block
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
                    # logic for finding the index of the ref positions
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


def find_all_amino_acids_fast(file_path):
    amino_acids_block = {}
    amino_acids_counter = []
    seq_count = 0
    seq_index = 0
    position_index = 0
    position_indices = []
    line_start = 0
    reference = []
    index_dict = []
    with open(file_path, 'r') as file:
        for line in file:
            # get line start information as well as the index of the position
            if line.__contains__('Prot'):
                if len(amino_acids_block) > 0:
                    for i, aa in amino_acids_block.items():
                        amino_acids_counter.append(Counter(aa))
                    amino_acids_block = {}
                line_start = line.index(line.split()[1][0])
                line_start_index = int(line.split()[1])
                seq_index = 0
                seq_count = seq_index
                next(file)
            elif re.search('[A-Za-z]+\*[0-9]+:[0-9]+.*', line):
                # add reference at start of block
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

                    # logic for finding the index of the ref positions
                    position_indices = []
                    position_index = line_start
                    i = 0
                    ## with or without +1??? without more variation but counting looks like +1 is valid??
                    while i < len(ref):
                        a = line[position_index:]
                        if not (line[position_index] == ' ' or line[position_index] == '.'):
                            i += 1
                            position_indices.append(position_index)
                            amino_acids_block[position_index] = [line[position_index]]
                        position_index += 1
                else:
                    for i in position_indices:
                        if i < len(line):
                            if line[i] == '-':
                                amino_acids_block[i].append(amino_acids_block[i][0])
                            elif line[i] == '\n':
                                amino_acids_block[i].append('.')
                                # print(len(line))
                            elif line[i] == ' ':
                                amino_acids_block[i].append('.')
                                # print(len(line))
                            else:
                                amino_acids_block[i].append(line[i])
                        else:
                            amino_acids_block[i].append('.')
                seq_index += 1

    if len(amino_acids_block) > 0:
        for i, aa in amino_acids_block.items():
            amino_acids_counter.append(Counter(aa))

    return amino_acids_counter
