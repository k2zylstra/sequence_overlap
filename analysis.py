'''
Uses the data generated by the compression_data.py file to generate graphs
and perform analysis
'''

import math
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

PATH = "/home/kieran/Documents/Bachelor_Thesis/Implementation/" # working directory path
GENE_FILE = "hg38.ensGene.gtf" # gene data file
RESULTS_FILE = "compression_rates_results.csv" # main data file
RESULTS_DELIM_FILE = "compression_rates_results_delim.csv" # output delimiter file
LABELS = ["Control", "Static Reduction", "Static Reduction with Seperate Lists", \
    "Dynamic Reduction with Deliminator", "Augmented with Index", "Augmented with Small Header Table"] # method labels

# creates the histogram graph of all methods listed in M.
def create_graph_hist(file_path,file_name, labels):
    full_path = file_path+file_name
    df = pd.read_csv(full_path, index_col=0)

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

# creates the pdf graph of all methods listed in M.
def create_graph_hist_pdf(file_path,file_name, labels):
    full_path = file_path+file_name
    df = pd.read_csv(full_path, index_col=0)

    # methods to compare
    M = [2,3,4,5]
    colors = ["blue", "coral", "green", "purple", "red", "yellow"]
    df_clean = pd.DataFrame()

    fig, ax  = plt.subplots(figsize=(15,5))

    for m in zip(M,colors):
        df0 = df[df["method_num"] == m[0]]
        for s in df["set_size"].unique():
            df_set = df0[df0["set_size"] == s]
            df_set["space"] = df_set["space"]/s
            df_clean = df_clean.append(df_set, ignore_index=True)
        df_clean[df_clean["method_num"] == m[0]]["space"].plot.hist(ax=ax, bins=20,alpha=.7,color=m[1], density=True)
        df_clean[df_clean["method_num"] == m[0]]["space"].plot.kde(ax=ax, color=m[1])
    df_clean.to_csv("results_normalized_histogram.csv")
    
    handles = [Rectangle((0,0),1,1,color=c,ec="k") for c in colors]

    # selects wich labels to be used
    l = []
    for m in M:
        l.append(labels[m])

    plt.legend(handles, l, loc="upper left")
    ax.set_ylabel("Probability per Bit")
    ax.set_xlabel("Bits Taken on Average per Number")
    plt.title("Histogram of Space Taken for Each Method")
    plt.show()
    return

# this creates the trendline graph of all methods listed in M
#   as the set size increases
def create_graph_trendline(file_path, file_name, labels):
    full_path = file_path+file_name
    df = pd.read_csv(full_path,index_col=0)

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
    

# This creates the bar graph that lists the average space
#   taken by all methods listed in M
def create_graph_baravg(file_path, file_name, labels):
    full_path = file_path+file_name
    df = pd.read_csv(full_path,index_col=0)
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

# this creates the bar graph for method 2 demonstrating how the 
#   method scales with different set sizes
def create_graph_bar_meth2(file_path, file_name, labels):
    full_path = file_path + file_name
    df = pd.read_csv(full_path, index_col=0)
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
       
# creates the line graph for methods 3, 4, and 5 to showcase how they change
#   based on dilimter size
def create_graph_lin_delim(file_path, file_name):
    full_path = file_path + file_name
    df = pd.read_csv(full_path, index_col=0)

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

def main():
    create_graph_hist(PATH, RESULTS_FILE, LABELS)
    create_graph_hist_pdf(PATH, RESULTS_FILE, LABELS)
    create_graph_trendline(PATH, RESULTS_FILE, LABELS)
    create_graph_baravg(PATH, RESULTS_FILE, LABELS)
    create_graph_bar_meth2(PATH, RESULTS_FILE, LABELS)
    create_graph_lin_delim(PATH, RESULTS_DELIM_FILE)

if __name__ == "__main__":
    main()