'''
Downloads genomic features and then runs compression tests
'''

import requests, json
import math
import pandas as pd
import numpy
import random
import check_prob_bitstring
import matplotlib.pyplot as plt
import seaborn as sns

# will contain repeats
def create_rnd_rows(num_rows, file_len):
    rndlist = []
    for i in range(num_rows):
        rndlist.append(random.random() * file_len // 1)
    
    rndlist.sort()
    return rndlist


def get_data_ucsc_gtf(num_rows, file_name):
    #create empty data frame with gtf format
    df = pd.DataFrame(columns=['sequence', 'source', 'feature', 'start', 'end'])
    #get file length
    num_lines = sum(1 for line in open(file_name, 'r'))

    rows = create_rnd_rows(num_rows, num_lines)
    #read in specified rows
    with open(file_name, 'r') as f:
        df_count = 0
        for i, line in enumerate(f):
            if i in rows:

                df.loc[df_count] = line.split()[:5]
                df_count += 1
    df = df.drop(columns=['sequence', 'source', 'feature'])
    df = df.astype(int)
    arr = df.to_numpy()
    return arr

def get_data_encode_api():
    headers = {'accept': 'application/json'}

    search_terms = ['bone', 'chip']
    search_phrase = ""
    search_terms_s = len(search_terms)
    for i in range(search_terms_s):
        search_phrase += search_terms[i]
        if i == search_terms_s - 1:
            break
        search_phrase += "+"

    url = 'https://www.encodeproject.org/search/?searchTerm=' + search_phrase + '&frame=object'
    print(url)

    response = requests.get(url, headers=headers)
    
    results = response.json()

    print(json.dumps(results, indent=4))

# don't know how to put my methodology into computer form yet
#   so I am putting the values in manually temporarily
def get_prob_nat_delim(delim_size, avg_num_size):
    p = check_prob_bitstring.get_prob(avg_num_size, delim_size)
    return p


# takes in double array of start and end values
def test_comp_normal(A):

    # all this does is figure 32 bits for each number
    #   this is more a formality as a function than real use
    return len(A) * 32


# this is the version that has just a smaller base int size used
# takes in double array of start and end values
def test_comp_first(A):

    # finds the smallest size that it can store the bits in
    largest_bit_size = 0
    for i in range(len(A)):
        for j in range(len(A[i])):
            bits = math.ceil(math.log2(A[i][j]))
            if bits > largest_bit_size:
                largest_bit_size = bits
    return largest_bit_size * len(A)

# this splits all the numbers to seperate lists so the maximum
#   universe size is the maximum length of a gene. In practice
#   this wil be larger because this is assuming that every number
#   gets its own list which is not practical, but it is a lower
#   bound
def test_comp_second(A):
    largest_bit_size = 0
    for i in range(len(A)):

        # this is start - end
        num = A[i][1] - A[i][0] 
        if num == 0:
            continue
        bits = math.ceil(math.log2(num))
        if bits > largest_bit_size:
            largest_bit_size = bits
    return largest_bit_size * len(A)

# this uses a delimiter of size 'delim_size' to separate all items
#   in a continual bit-string. This also assumes the seperate lists
#   as in the above function
def test_comp_delim(A, delim_size):
    size = 0

    avg = 0
    for i in range(len(A)):

        # this is start - end
        num = A[i][1] - A[i][0] 
        if num == 0:
            continue
        avg += num
        size += math.ceil(math.log2(num))
        size += delim_size
    
    avg = avg / len(A)
    avg_num_size = math.ceil(math.log2(avg))
    
    # get probability that delim naturally occurs
    prob_nat_delim = get_prob_nat_delim(delim_size, avg_num_size) 
    size += prob_nat_delim * len(A) * 32
    return size

# this uses two delimiters to augment the data structure and also
#   store the index of the item
def test_comp_delim2(A, delim_size):
    size = 0
    avg = 0

    for i in range(len(A)):

        # this is start - end
        num = A[i][1] - A[i][0]
        if num == 0:
            continue
        avg += num
        size += math.ceil(math.log2(num))
        size += 2 * delim_size
    
    avg = avg / len(A)
    avg_num_size = math.ceil(math.log2(avg))
    
    # gets probability that delim naturally occurs
    prob_nat_delim = get_prob_nat_delim(delim_size, avg_num_size) 
    
    # the two comes from the fact that the delim is twice as likely to occur
    size += 2 * prob_nat_delim * len(A) * 32
    return size

# this uses the two delimiters method, but also tries to minimize
#   the header table storing all naturally occuring delims by
#   making all spaces used only as long as the longest int
def test_comp_delim2_smlhdr(A, delim_size):
    size = 0
    avg = 0

    for i in range(len(A)):

        # this is start - end
        num = A[i][1] - A[i][0]
        if num == 0:
            continue
        avg += num
        size += math.ceil(math.log2(num))
        size += 2 * delim_size
    
    avg = avg / len(A)
    avg_num_size = math.ceil(math.log2(avg))
    
    # get probability that delim naturally occurs
    prob_nat_delim = get_prob_nat_delim(delim_size, avg_num_size)

    # this assumes that the largest index is what gets stored
    #   in the header table. That is not always true. To find
    #   out we would need to calculate which numbers actually
    #   contained the delimiters and use their indexes, but this
    #   is a good estimate
    largest_bit_size = math.ceil(math.log2(len(A) - 1))
    size += 2 * prob_nat_delim * len(A) * largest_bit_size
    return size


def main():

    # for sequence set sizes of variable, 'set_size'
    set_size = [100, 1000, 10000]

    # holds data for all found sizes
    S = {}
    
    for s in set_size:
        size = []
        A = get_data_ucsc_gtf(s, "/home/kieran/Documents/Bachelor_Thesis/Implementation/hg38.ensGene.gtf")
        size.append(test_comp_normal(A))
        size.append(test_comp_second(A))
        size.append(test_comp_delim(A, 9))
        size.append(test_comp_delim2(A, 9))
        size.append(test_comp_delim2_smlhdr(A, 9))

        S[s] = size
    
    with open("compression_rates_results.csv", "w") as f:
        for i in range(len(set_size)):
            f.write(str(set_size[i]))
            if i == len(set_size) - 1:
                break
            f.write(",")
        f.write("\n")
        for i in set_size:
            for j in range(len(S[i])):
                f.write(str(S[i][j]))
                if i == len(set_size) - 1:
                    break
                f.write(",")
            f.write("\n")

    plt.bar(x=S[100], height=5, width=0.8)
    plt.show()
    
def create_graph(file_path):
    sns.set_style("dark")
    df = pd.read_csv(file_path)
    df = df.set_index('set_size')
    df.plot.bar()
    plt.show()
    return

if __name__ == "__main__":
    #main()
    create_graph("/home/kieran/Documents/Bachelor_Thesis/Implementation/compression_rates_results.csv")