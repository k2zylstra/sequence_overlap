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
from matplotlib.patches import Rectangle
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
     
def create_graph_hist(file_path, labels):
    df = pd.read_csv(file_path, index_col=0)

    # methods to compare
    M = [0,1,2,3,4,5]
    colors = ["blue", "coral", "green", "purple", "red", "yellow"]
    df_clean = pd.DataFrame()

    fig, ax  = plt.subplots(figsize=(15,5))
    #fig, ax  = plt.subplots(figsize=(9,5))

    for m in zip(M,colors):
        df0 = df[df["method_num"] == m[0]]
        for s in df["set_size"].unique():
            df_set = df0[df0["set_size"] == s]
            df_set["space"] = df_set["space"]/s
            df_clean = df_clean.append(df_set, ignore_index=True)
        df_clean[df_clean["method_num"] == m[0]]["space"].plot.hist(ax=ax, bins=20,alpha=.7,color=m[1], density=False)
        #df_clean[df_clean["method_num"] == m[0]]["space"].plot.kde(ax=ax, color=m[1])
    df_clean.to_csv("results_normalized_histogram.csv")
    
    handles = [Rectangle((0,0),1,1,color=c,ec="k") for c in colors]

    # selects wich labels to be used
    l = []
    for m in M:
        l.append(labels[m])

    plt.legend(handles, l, loc="upper left")
    ax.set_ylabel("Count")
    #ax.set_ylabel("Probability per Bit")
    ax.set_xlabel("Bits Taken on Average per Number")
    plt.ylim(0, 1700)
    plt.title("Histogram of Space Taken for Each Method")
    plt.show()
    return

def create_graph_trendline(file_path, labels):
    df = pd.read_csv(file_path,index_col=0)

    # take average
    # condense into 6 x 3 points (18)
    # plot

    # M stands for method. These are the methods that will be plotted.
    #   Change to compare individual methods
    M = [0,1,2,3,4,5]

    df_clean = pd.DataFrame()
    for i in df["set_size"].unique():
        for m in M:
            df0 = df.loc[df["set_size"] == i]
            df0 = df0.loc[df0["method_num"] == m]
            df0["space"] = df0["space"].mean()
            df0 = df0.drop_duplicates()
            df_clean =df_clean.append(df0)
    
    fig, ax  = plt.subplots(figsize=(9,5))

    for m in M:
        df_clean.loc[df_clean["method_num"] == m].plot(kind="line", x="set_size", y="space", ax=ax, label=labels[m])
        print(df_clean.loc[df_clean["method_num"] == m].head())

    ax.set_ylabel("Bits Taken")
    ax.set_xlabel("Set Size")
    plt.title("Trendline of Increasing Set Size")

    plt.show()
    

def create_graph_baravg(file_path, labels):
    df = pd.read_csv(file_path,index_col=0)
    M = [0,1,2,3,4,5]
    columns = ["method_num", "avg"]
    
    df_norm = pd.DataFrame()
    df_clean = pd.DataFrame()
    l = []
    for m in M:
        df0 = df.loc[df["method_num"] == m]
        for s in df["set_size"].unique():
            df_set = df0.loc[df0["set_size"] == s]
            df_set["space"] = df_set["space"]/s
            df_set.drop(["trial_num", "set_size"], axis=1, inplace=True)
            df_set = df_set.drop_duplicates()
            df_norm = df_norm.append(df_set, ignore_index=True)
    for m in M:
        df0 = df_norm.loc[df_norm["method_num"] == m]
        df0["space"] = df0["space"].mean()
        df0 = df0.drop_duplicates()
        df_clean = df_clean.append(df0, ignore_index=True)

        l.append(labels[m])
    
    #df_clean = df_clean.assign(names=l)
    fig, ax = plt.subplots(figsize=(9,5))
    df_clean.plot.bar(x="method_num", y="space", ax=ax)    
    ax.set_xlabel("Method")
    ax.set_ylabel("Average Bits Used per Number")
    plt.title("Average Bits Used")
    x = df_clean["method_num"].tolist()
    y = df_clean["space"].tolist()
    for i in range(len(y)):
        y[i] = math.trunc(y[i] * 1000) / 1000
    for i in range(len(x)):
        plt.text(i, y[i], y[i], ha="center")
    plt.show()

def create_graph_bar_meth2(file_path, labels):
    df = pd.read_csv(file_path, index_col=0)
    df = df.loc[df["method_num"] == 2]
    df = df.drop("method_num",axis=1)
    df = df.drop("trial_num",axis=1)
    for s in df["set_size"].unique():
        df.loc[df["set_size"] == s, "space"] = df.loc[df["set_size"] == s]["space"].mean()
    df = df.drop_duplicates()
    for s in df["set_size"].unique():
        df.loc[df["set_size"] == s, "space"] = df.loc[df["set_size"] == s]["space"].div(s)
    
    fig, ax = plt.subplots(figsize=(9,6))
    df.plot.bar(x="set_size", y="space", legend=False, ax=ax)
    
    ax.set_xlabel("Set Size")
    ax.set_ylabel("Average Bits Used per Number")
    plt.title("Average Bits Used in Method 2 Depending on Set Size")

    x = df["set_size"].tolist()
    y = df["space"].tolist()
    for i in range(len(y)):
        y[i] = math.trunc(y[i] * 1000) / 1000
    for i in range(len(x)):
        plt.text(i, y[i], y[i], ha="center")

    plt.show()

def mean_number_values(file_path):
    A = get_data_ucsc_gtf(100000, file_path)
    total = 0
    for a in A:
        total += a[1] - a[0]
    avge = total / len(A)
    print("The mean length of genomic feature is:",avge)
        
def graph_scatter_delim(file_path):
    df = pd.read_csv(file_path, index_col=0)

    for d in df["delim_size"].unique():
        for m in df["method_num"].unique():
            df.loc[df["method_num"].eq(m) & df["delim_size"].eq(d), "space"] = df.loc[df["method_num"].eq(m) & df["delim_size"].eq(d), "space"].mean()
    df = df.drop("trial_num", axis=1)
    df = df.drop_duplicates()

    fig, ax = plt.subplots(figsize=(9,5))
    colors = ["red", "blue", "green"]
    labels = ["Dynamic Reduction with Deliminator", "Augmented with Index", "Augmented with Small Header Table"]
    for m in zip(df["method_num"].unique(), colors, labels):
        df.loc[df["method_num"] == m[0]].plot(kind="line", x="delim_size", y="space", c=m[1], ax=ax, label=m[2])
        print(df.loc[df["method_num"] == m[0]])
    ax.set_ylabel("Space Used (bits)")
    ax.set_xlabel("Deiliminator Size (bits)")
    plt.title("Space Used by Deliminator Size")
    
    plt.show()
#======================================================================== 

def test_delim_size(results_file):
    PATH = "/home/kieran/Documents/Bachelor_Thesis/"
    GENE_FILE = "hg38.ensGene.gtf"
    RESULTS_FILE = "compression_rates_results_delim.csv"
    num_trials = 300 
    set_size = 500
    delim_range = range(0,16)

    now = datetime.datetime.now()
    print(now)

    columns = ["trial_num", "method_num", "delim_size", "space"]
    df = pd.DataFrame(columns=columns)
    for d in delim_range:
        for i in range(num_trials):
            A = get_data_ucsc_gtf(set_size, PATH+GENE_FILE)
            trial = []
            trial.append(test_comp_delim(A, d))
            trial.append(test_comp_delim2(A, d))
            trial.append(test_comp_delim2_smlhdr(A, d))
            for t in range(len(trial)):
                row = pd.Series([i, t, d, trial[t]],index=df.columns)
                df = df.append(row, ignore_index=True)
        now = datetime.datetime.now()
        print("Finished trials for delimiter size:", str(d), ": at", str(now))
    
    print(df.head())
    df.to_csv(RESULTS_FILE)
    
    now = datetime.datetime.now()
    print(now)



def main():
    PATH = "/home/kieran/Documents/Bachelor_Thesis/"
    GENE_FILE = "hg38.ensGene.gtf"
    RESULTS_FILE = "compression_rates_results.csv"
    num_trials = 1000 

    now = datetime.datetime.now()
    print(str(now))

    # for sequence set sizes of variable, 'set_size'
    set_sizel = [100, 1000, 10000]

    # holds data for all found sizes
    # Is a dictionary with keys as the set sizes and values as a 2d array
    #   of columns being the methods and rows being each trial
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
    file_path = "/home/kieran/Documents/Bachelor_Thesis/Implementation/compression_rates_results.csv"
    file_path = "compression_rates_results.csv"
    labels = ["Control", "Static Reduction", "Static Reduction with Seperate Lists", "Dynamic Reduction with Deliminator", "Augmented with Index", "Augmented with Small Header Table"]
    create_graph_hist("compression_rates_results.csv", labels)
    #main()
    #create_graph_trendline(file_path, labels)
    #create_graph_baravg(file_path, labels)
    #test_delim_size("delim_results.csv")
    create_graph_bar_meth2(file_path, labels)


    PATH = "/home/kieran/Documents/Bachelor_Thesis/"
    GENE_FILE = "hg38.ensGene.gtf"
    #mean_number_values(PATH+GENE_FILE)
    #graph_scatter_delim("compression_rates_results_delim.csv")
