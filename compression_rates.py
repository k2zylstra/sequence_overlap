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
import os
import datetime

# will contain repeats
def create_rnd_rows(num_rows, file_len):
    rndlist = []
    i = 0
    while i < num_rows:
        rnd = random.random() * file_len // 1
        if rnd in rndlist:
            continue
        rndlist.append(rnd)
        i+=1
    
    rndlist.sort()
    return rndlist


def get_data_ucsc_gtf(num_rows, file_name):
    #create empty data frame with gtf format
    df = pd.DataFrame(columns=['sequence', 'source', 'feature', 'start', 'end'])
    #get file length
    #num_lines = sum(1 for line in open(file_name, 'r'))
    num_lines = 0
    with open(file_name, 'r') as f:
        for i, line in enumerate(f):
            num_lines += 1

    rows = create_rnd_rows(num_rows, num_lines)
    num_rows = len(rows)
    #read in specified rows
    with open(file_name, 'r') as f:
        df_count = 0
        for i, line in enumerate(f):
            if i == rows[0]:
                rows.pop(0)
                df.loc[df_count] = line.split()[:5]
                df_count += 1
                if not rows:
                    break
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

def create_new_csv(path, num_trials):
    if os.path.exists(path):
        with open("compression_rates_results.csv", "w") as f:
            f.write("set_size|")
            for i in range(num_trials):
                f.write("Trial"+str(i+1))
                if i == num_trials-1:
                    break
                f.write("|")
            f.write("\n")

# this takes in the dictionary filled with set size as the key and
#   the array filled with the results of each method as the value. From
#   this it transforms the data to a stacked verion of this with added columns
#   of trial number and method number and returns a data frame
def transform_data(data):
    new_df = pd.DataFrame(columns=["set_size", "trial_num", "method_num", "space"])
    row = []

    # for each dictionary key in data
    for dk in data:
        for i in range(len(data[dk])):
            for j in range(len(data[dk][i])):
                row.append(dk)
                row.append(i)
                row.append(j)
                row.append(data[dk][i][j])
                s_row = pd.Series(row, index=new_df.columns)
            
                new_df = new_df.append(s_row,ignore_index=True)
                row = []
    
    return new_df



# Data is a dictionary with set sizes as keys and a 2d array of compression type as rows
#   and the amount of space as the columns
def write_csv_old(path, resultsfile, data, set_sizel, num_trials):

    if not os.path.exists(path+resultsfile):
        create_new_csv(path, num_trials)
    with open("compression_rates_results.csv", "a") as f:
        for i in set_sizel:
            f.write(str(i)+"|")
            for j in range(len(data[i])):
                f.write(str(data[i][j]))
                if j == len(data[i]) - 1:
                    break
                f.write("|")
            f.write("\n")
     
def create_graph(file_path):
    df = pd.read_csv(file_path, index_col=0)

    # methods to compare
    M = [1]
    colors = ["blue", "coral", "green", "purple", "red"]
    df_clean = pd.DataFrame()

    for m in zip(M,colors):
        df0 = df[df["method_num"] == m[0]]
        for s in df["set_size"].unique():
            df_set = df0[df0["set_size"] == s]
            df_set["space"] = df_set["space"]/s
            df_clean = df_clean.append(df_set, ignore_index=True)
        df_clean[df_clean["method_num"] == m[0]]["space"].plot.hist(bins=20,alpha=.5,color=m[1])

    plt.show()
    return

#======================================================================== 

def main():
    PATH = "/home/kieran/Documents/Bachelor_Thesis/Implementation/"
    GENE_FILE = "hg38.ensGene.gtf"
    RESULTS_FILE = "compression_rates_results.csv"
    num_trials = 1000 

    now = datetime.datetime.now()
    print(str(now))

    # for sequence set sizes of variable, 'set_size'
    set_sizel = [100, 1000, 10000]

    # holds data for all found sizes
    S = {}
    
    for s in set_sizel:
        size = []
        for i in range(num_trials):
            A = get_data_ucsc_gtf(s, PATH+GENE_FILE)

            trial = []
            trial.append(test_comp_normal(A))
            trial.append(test_comp_first(A))
            trial.append(test_comp_second(A))
            trial.append(test_comp_delim(A, 9))
            trial.append(test_comp_delim2(A, 9))
            trial.append(test_comp_delim2_smlhdr(A, 9))
            
            #now = datetime.datetime.now()
            #print("completed trial:", str(i),"| timestamp:", str(now))
            if i % 100 == 0:
                print("Completed trial", str(i), "for set size", str(s))

            size.append(trial)
        S[s] = size
        now = datetime.datetime.now()
        print("completed trials for set size:", str(s), "| timestamp:", str(now))

    df = transform_data(S)
    print(df.head())
    df.to_csv(PATH+RESULTS_FILE)
    #write_csv_old(PATH, RESULTS_FILE, S, set_sizel, num_trials)
    now = datetime.datetime.now()
    print(str(now))
    #create_graph(PATH+RESULTS_FILE)
 
if __name__ == "__main__":
    main()
    #create_graph("/home/kieran/Documents/Bachelor_Thesis/Implementation/compression_rates_results.csv")