import sys
import time
### change this to your path ###
sys.path.append('/Users/danielwang/Desktop/CM122/cm122_project1b')

def load_reference(file_path):
    with open(file_path) as f:
        reference_genome = f.read().splitlines()
    reference_genome.pop(0)
    return "".join(reference_genome)

def load_single_reads(file_path):
    with open(file_path) as f:
        reads = f.read().splitlines()
    return [read for read in reads if read[0] != '>']

def load_paired_reads(file_path):
    with open(file_path) as f:
        reads = f.read().splitlines()
    return [[reads[i], reads[i+2]] for i in range(1, len(reads), 4)]

class suffix:
     
    def __init__(self):
         
        self.index = 0
        self.rank = [0, 0]

def suffix_array(genome):
    n = len(genome)
    suffixes = [suffix() for i in range(n)]
 
    for i in range(n):
        suffixes[i].index = i
        suffixes[i].rank[0] = (ord(genome[i]) -
                               ord("a"))
        suffixes[i].rank[1] = (ord(genome[i + 1]) -
                        ord("a")) if ((i + 1) < n) else -1
 
    suffixes = sorted(
        suffixes, key = lambda x: (
            x.rank[0], x.rank[1]))
    ind = [0] * n
    k = 4
    while (k < 2 * n):
        rank = 0
        prev_rank = suffixes[0].rank[0]
        suffixes[0].rank[0] = rank
        ind[suffixes[0].index] = 0
 
        for i in range(1, n):
            if (suffixes[i].rank[0] == prev_rank and
                suffixes[i].rank[1] == suffixes[i - 1].rank[1]):
                prev_rank = suffixes[i].rank[0]
                suffixes[i].rank[0] = rank
            else: 
                prev_rank = suffixes[i].rank[0]
                rank += 1
                suffixes[i].rank[0] = rank
            ind[suffixes[i].index] = i
 
        for i in range(n):
            nextindex = suffixes[i].index + k // 2
            suffixes[i].rank[1] = suffixes[ind[nextindex]].rank[0] \
                if (nextindex < n) else -1
 
        suffixes = sorted(
            suffixes, key = lambda x: (
                x.rank[0], x.rank[1]))
        k *= 2
    suffixArr = [0] * n
     
    for i in range(n):
        suffixArr[i] = suffixes[i].index
    return suffixArr

def bwt_from_suffix(genome, s_array=None):
    if s_array is None:
        s_array = suffix_array(genome)
    return([genome[idx - 1] for idx in s_array])

def bwt(genome):
    rotations = []
    length = len(genome)
    for i in range(length):
        rotations.append((str(genome[i:length] + genome[0:i]), i))
    rotations_sorted = sorted(rotations, key=lambda x: x[0])
    return [x[0][-1] for x in rotations_sorted], [x[1] for x in rotations_sorted]

def first_occurence(bwt_seq):
    sorted_bwt = sorted(bwt_seq)
    return {symbol: sorted_bwt.index(symbol) for symbol in ['A', 'G', 'C', 'T']}

def count_symbols(bwt_seq):
    counts = {'A': [0], 'G': [0], 'C': [0], 'T': [0]}
    for i in range(len(bwt_seq)):
        symbol = bwt_seq[i]
        for base in counts:
            if base == symbol:
                counts[base].append(counts[base][i] + 1)
            else:
                counts[base].append(counts[base][i])
    return counts


def bw_matching(bwt_seq, patterns, suf_arr, first_occs, symbol_counts):
    match_locations = {}
    for pattern in patterns:
        top = 0
        bottom = len(bwt_seq) - 1
        pattern_og = pattern
        while top <= bottom:
            if pattern:
                symbol = pattern[-1]
                pattern = pattern[0:len(pattern)-1]
                positions = []
                for i in range(top, bottom+1):
                    if bwt_seq[i] == symbol:
                        positions.append(i)
                        break
                if positions:
                    top = first_occs[symbol] + symbol_counts[symbol][top]
                    bottom = first_occs[symbol] + symbol_counts[symbol][bottom+1] - 1
                else:
                    match_locations[pattern_og] = []
                    break
            else:
                match_pos = suf_arr[top:bottom+1]
                match_locations[pattern_og] = match_pos
                break
    return match_locations

def split_read(read, k, l):
    i = 0
    k_mers = []
    while i + k <= l:
        k_mers.append(read[i:i+k])
        i += k
    if i != l:
        leftover_len = l - i
        if leftover_len < 6:
            k_mers[-1] += read[i:]
        else:
            k_mers.append(read[i:])
    return k_mers

class Mutation:
    def __init__(self, type, position, base_1, base_2=''):
        self.type = type
        self.position = position
        if not base_2:
            self.description = '>' + type + str(position) + ' ' + base_1
        else:
            self.description = '>' + type + str(position) + ' ' + base_1 + ' ' + base_2
    
    def __hash__(self):
        return hash((self.type, self.position, self.description))

    def __eq__(self, other):
        if not isinstance(other, type(self)): return NotImplemented
        return self.type == other.type and self.position == other.position and self.description == other.description

def find_substitution(reference_seq, read_seq, start_idx):
    mutations = []
    max_subs = 1
    num_subs = 0
    for i in range(len(reference_seq)):
        b_ref = reference_seq[i]
        b_read = read_seq[i]
        if b_ref != b_read:
            num_subs += 1
            if num_subs <= max_subs:
                
                mutations.append(Mutation('S', start_idx + i, b_ref, b_read))
            else:
                mutations.clear()
                break
    return mutations

def find_indel(reference_seq, read_seq, start_idx):
    mutations = []
    i = 0
    ref_len = len(reference_seq)
    read_len = len(read_seq)
    shorter_seq_len = min([ref_len, read_len])
    while i < shorter_seq_len:
        b_ref = reference_seq[i]
        b_read = read_seq[i]
        if b_ref != b_read:
            if ref_len < read_len: # check for insertion
                read_new = read_seq[0:i]
                if i < read_len - 1:
                    read_new += read_seq[i+1:]
                if read_new == reference_seq:
                    mutations.append(Mutation('I', start_idx + i - 1, b_read))
                    break
            else: # check for deletion
                read_new = read_seq[0:i] + b_ref + read_seq[i:]
                if read_new == reference_seq:
                    mutations.append(Mutation('D', start_idx + i, b_ref))
                    break
        i += 1
    if not mutations:
        if ref_len < read_len: # insertion
            read_new = read_seq[0:i]
            if read_new == reference_seq:
                mutations.append(Mutation('I', start_idx + i - 1, read_seq[i]))
        else: # deletion
            read_new = read_seq + reference_seq[i]
            if read_new == reference_seq:
                mutations.append(Mutation('D', start_idx + i, reference_seq[i]))
    return mutations
                
reference = load_reference('project1b_1000000-data/project1b_1000000_reference_genome.fasta')

#bw_seq, suf_arr = bwt(reference)
reference += '$'
start = time.time()
suf_arr = suffix_array(reference)
end = time.time()
print('suffix array time: ' + str(end - start), file=sys.stderr)

start = time.time()
bw_seq = bwt_from_suffix(reference, suf_arr)
first_occs = first_occurence(bw_seq)
symbol_counts = count_symbols(bw_seq)
end = time.time()
print('bwt time: ' + str(end - start), file=sys.stderr)

reads_path = 'project1b_1000000-data/project1b_1000000_with_error_paired_reads.fasta'

mutations = []
lines_read = 0
want_print = False
start = time.time()
with open(reads_path) as f:
    pair = []
    for line in f:
        if lines_read % 2 == 1:
            pair.append(line[0:len(line)-1])
        if lines_read % 4 == 3:
            prev_read_end = -1 # position where the first read in the pair ends after alignment
            for i in range(len(pair)):
                seq_read = pair[i]
                l = len(seq_read)
                k = 10
                k_mers = []
                if l > k:
                    k_mers = split_read(seq_read, k, l) # only last fragment will ever not have length k
                else:
                    k_mers.append(seq_read)
                n = len(k_mers)
                match_positions = bw_matching(bw_seq, k_mers, suf_arr, first_occs, symbol_counts)
                
                # preprocess matches and determine potential starting indices
                expected_starts = [] # predicted start of the full k-mer
                for j in range(n):
                    k_mer = k_mers[j]
                    valid_positions = []
                    for matched_pos in match_positions[k_mer]:
                        if (matched_pos > prev_read_end and i == 1) or i == 0:
                            expected_starts.append(matched_pos - j * k)
                            valid_positions.append(matched_pos)
                    match_positions[k_mer] = valid_positions

                for st in expected_starts:
                    expected = st
                    unaligned_k_mers = []
                    potential_mutations = []
                    for x in range(n):
                        k_mer = k_mers[x]
                        k = len(k_mer)
                        affected_by_mutation = False
                        if match_positions[k_mer]:
                            correctly_sequenced = False
                            for mp in match_positions[k_mer]:
                                if mp == expected:
                                    expected += k
                                    correctly_sequenced = True
                                    break
                            if not correctly_sequenced:
                                affected_by_mutation = True
                        else:
                            affected_by_mutation = True
                        
                        if affected_by_mutation:
                            if want_print:
                                print("current k_mer: " + k_mer, file=sys.stderr)
                                print("expected: " + str(expected), file=sys.stderr)
                            
                            new_mutation = find_substitution(reference[expected:expected+k], k_mer, expected)
                            if new_mutation:
                                
                                if want_print:
                                    print([x.description for x in new_mutation], file=sys.stderr)
                                    print(match_positions, file=sys.stderr)
                                    print(reference[expected:expected+k], file=sys.stderr)
                                    print('\n', file=sys.stderr)
                                
                                potential_mutations += new_mutation
                                expected += k
                            else:
                                if x == 0:
                                    insertion_start = expected + 1
                                    deletion_start = expected - 1

                                    new_mutation = find_indel(reference[insertion_start:expected+k], k_mer, insertion_start) # insertion
                                    if new_mutation:
                                        if want_print:
                                            print([x.description for x in new_mutation], file=sys.stderr)
                                            print(match_positions, file=sys.stderr)
                                            print(reference[insertion_start:expected+k], file=sys.stderr)
                                            print('\n', file=sys.stderr)
                                        
                                        potential_mutations += new_mutation
                                    else:
                                        if deletion_start > 0:
                                            new_mutation = find_indel(reference[deletion_start:expected+k], k_mer, deletion_start) # deletion
                                        else:
                                            new_mutation = []
                                        if new_mutation:
                                            if want_print:
                                                print([x.description for x in new_mutation], file=sys.stderr)
                                                print(match_positions, file=sys.stderr)
                                                print(reference[deletion_start:expected+k], file=sys.stderr)
                                                print('\n', file=sys.stderr)
                                            
                                            potential_mutations += new_mutation
                                    expected += k
                                else:
                                    new_mutation = find_indel(reference[expected:expected+k-1], k_mer, expected) # insertion
                                    if new_mutation:
                                        
                                        if want_print:
                                            print([x.description for x in new_mutation], file=sys.stderr)
                                            print(match_positions, file=sys.stderr)
                                            print(reference[expected:expected+k-1], file=sys.stderr)
                                            print('\n', file=sys.stderr)
                                        
                                        potential_mutations += new_mutation
                                        expected += k - 1
                                    else:
                                        if expected + k + 1 < len(reference):
                                            new_mutation = find_indel(reference[expected:expected+k+1], k_mer, expected) # deletion
                                        else:
                                            new_mutation = []
                                        if new_mutation:
                                            
                                            if want_print:
                                                print([x.description for x in new_mutation], file=sys.stderr)
                                                print(match_positions, file=sys.stderr)
                                                print(reference[expected:expected+k+1], file=sys.stderr)
                                                print('\n', file=sys.stderr)
                                            
                                            potential_mutations += new_mutation
                                            expected += k + 1
                                        else:
                                            
                                            if want_print:
                                                print("unsuccessful alignment of: " + k_mer, file=sys.stderr)
                                                print(match_positions, file=sys.stderr)
                                                if expected + k + 1 < len(reference):
                                                    print(reference[expected:expected+k+1], file=sys.stderr)
                                                else:
                                                    print(reference[expected:expected+k], file=sys.stderr)
                                                print('\n', file=sys.stderr)
                                            
                                            unaligned_k_mers.append(x)
                                            expected += k
                    if not unaligned_k_mers:
                        mutations += potential_mutations
                        prev_read_end = expected - 1
                        break
            pair = []
        lines_read += 1
end = time.time()
print('alignment time: ' + str(end-start), file=sys.stderr)

start = time.time()
# remove uncommon mutations
f_mutations = []
mutation_threshold = 2 # minimum number of occurrences of mutation
mutation_counts = {}
for mutation in mutations:
    if mutation not in mutation_counts:
        mutation_counts[mutation] = 1
    else:
        mutation_counts[mutation] += 1
mutations = set([m for m in mutation_counts if mutation_counts[m] >= mutation_threshold])

# remove substitutions near indels
mutation_types = [x.type for x in mutations]
mutation_positions = [x.position for x in mutations]
f_mutations = []
for z in range(len(mutations)):
    mutation = list(mutations)[z]
    if mutation_types[z] == 'S':
        pos = mutation_positions[z]
        count = 0
        for p in mutation_positions:
            if p in range(pos-2, pos+2):
                count += 1
        if count == 1:
            f_mutations.append(mutation)
    else:
        f_mutations.append(mutation)
end = time.time()
print('mutation adjustment time: ' + str(end - start), file=sys.stderr)

print("\n".join(sorted([x.description for x in f_mutations])))