import re
from collections import Counter


def find_amino_acids(file_path, position):
    """Finds amino acids and protein ids of a specific position in an alignment file"""
    amino_acids = []
    seq_index = 0
    position_index = 0
    line_start = 0
    reference = []
    index_dict = []
    amino_acids_ids = {}

    with open(file_path, 'r') as file:
        for line in file:
            # get line start information as well as the index of the position
            if line.__contains__('Prot'):
                line_start = line.index(line.split()[1][0])
                line_start_index = int(line.split()[1])
                seq_index = 0
                next(file)
            elif re.search('[A-Za-z]+\*[0-9]+:[0-9]+.*', line):
                # add reference at start of block
                if seq_index == 0:
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
                        while i < position - line_start_index:
                            if not (line[position_index] == ' ' or line[position_index] == '.'):
                                i += 1
                            position_index += 1

                        amino_acids.append(line[position_index])
                    else:
                        if position_index < len(line):
                            if line[position_index] == '-':
                                amino_acids.append(amino_acids[0])
                            elif line[position_index] == '\n':
                                amino_acids.append('.')
                            elif line[position_index] == ' ':
                                amino_acids.append('.')
                            else:
                                amino_acids.append(line[position_index])
                        else:
                            amino_acids.append('.')

                    prot_id = line.split()[0].split(':')[0]
                    if amino_acids_ids.__contains__(amino_acids[-1]):
                        if amino_acids_ids[amino_acids[-1]].__contains__(prot_id):
                            amino_acids_ids[amino_acids[-1]][prot_id] += 1
                        else:
                            amino_acids_ids[amino_acids[-1]][prot_id] = 1
                    else:
                        amino_acids_ids[amino_acids[-1]] = {prot_id: 1}

                elif position < line_start_index:
                    break
                seq_index += 1

    return amino_acids, amino_acids_ids


def find_all_amino_acids_fast(file_path):
    """Returns amino acids counts of all positions of the reference alignment in an alignment file"""
    amino_acids_block = {}
    amino_acids_counter = []
    seq_index = 0
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
                next(file)
            elif re.search('[A-Za-z]+\*[0-9]+:[0-9]+.*', line):
                # add reference at start of block
                if seq_index == 0:
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
                    while i < len(ref):
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


def cluster_counter(counter, aa_property, aa_annotations):
    """Clusters amino acids based on their annotations"""
    res = {}
    for aa, n in counter.items():
        if aa_annotations.__contains__(aa):
            if res.__contains__(aa_annotations[aa][aa_property]):
                res[aa_annotations[aa][aa_property]] += n
            else:
                res[aa_annotations[aa][aa_property]] = n
        else:
            res[aa] = n
    return res


def get_unique_proteins(amino_acids_ids, threshold=56, threshold_test=False):
    """Returns the unique protein id with a minimum amount higher than the set threshold,
    unless the threshold flag is set to True, then it will return the completeness and number of duplicates"""
    id_set = set()
    id_counted = []
    aa_ids_threshold = {}
    for aa, counts in amino_acids_ids.items():
        # print(aa)
        aa_ids_threshold[aa] = []
        for prot_id, c in counts.items():
            id_set.add(prot_id)
            if c > threshold:
                id_counted.append(prot_id)
                aa_ids_threshold[aa].append(prot_id)
                # print(prot_id)
    completeness = len(set(id_counted)) / len(id_set)
    duplicates = len(id_counted) - len(set(id_counted))
    if threshold_test:
        return completeness, duplicates
    else:
        return aa_ids_threshold


def calc_threshold_proteins(amino_acids_ids):
    """Function to analyze the amount of completeness and duplicates for different thresholds up to 200"""
    completeness_list = []
    duplicates_list = []
    for i in range(200):
        completeness, duplicates = get_unique_proteins(amino_acids_ids, i, True)
        completeness_list.append(completeness)
        duplicates_list.append(duplicates)
    return completeness_list, duplicates_list


def sort_by_other_dict(d, order):
    return {k: d[k] for k in order}
