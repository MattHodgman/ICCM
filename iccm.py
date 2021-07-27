'''

Process:

1. get cluster assignments df (already done in silhouette.py and celluster.py) -> find more clusters than optimal
2. essentially compute the jaccard indices for each pairwise comparison of clustering methods
    - thus for each pair of methods, compare all clusters against each other to find the number of cells that are in common divided by the number of total cells
    - jaccard index = intersection / union
    - then for each cluster, just pick the one in the other method that has the highest jaccard index
    - relabel the clusters for that method
3. now do the voting and assign cells to clusters if the majority of the methods agree on the cluster label
4. repeat on the points that are 'outliers' and then somehow combine them

'''


import argparse
import pandas as pd
from collections import Counter
import numpy as np

'''
Parse arguments.
'''
def parseArgs():
    parser = argparse.ArgumentParser(description='Implementation of the consensus clustering method proposed in https://www.sciencedirect.com/science/article/pii/S131915781930597X.')
    parser.add_argument('-i', '--input', help='Input cell cluster assignment files.', nargs='*', action='store', dest='input', required=True)
    parser.add_argument('-o', '--output', help='The directory to which output files will be saved.', type=str, required=False)
    parser.add_argument('-c', '--id', help='The name of the column that contains the item ID.', type=str, required=True)
    args = parser.parse_args()
    return args


'''
Get the header from an input csv file.
Returns a list where each element is a column header.
'''
def getHeader(file):
    # get first line of file aka the csv header
    with open(file) as f:
        header = f.readline().rstrip()
    f.close()

    # get the input file column headers as list
    header_columns = header.split(',')

    return header_columns


'''
Look at first line of an in input file (the header of the csv) to assess if it is the correct format.
'''
def validInput(file):

    # definition of the columns needed to constitute a valid input file
    NEEDED_COLUMNS = [id, CLUSTER, METHOD]

    input_columns = getHeader(file) # get input file header columns
    valid = all(col in input_columns for col in NEEDED_COLUMNS) # check if all needed columns are present in input columns list
    
    return valid


'''
Read file into dataframe.
'''
def getData(file):

    data = pd.read_csv(file, delimiter=',', index_col=id) # load data from input csv into dataframe
    method = data[METHOD].iloc[0] # get method name from first row (this is assuming at all rows are from the same method) NOTE: we may not want to assume this!!
    data = data.drop(METHOD, axis=1) # drop method column
    
    return data, method


'''
Compute the jaccard index of two clusters (sets of samples).
'''
def jaccard_index(cluster1, cluster2):
    i = len(cluster1.intersection(cluster2))
    u = len(cluster1.union(cluster2))
    jaccard_index = i / u
    return jaccard_index


'''
Get the method that has the fewest (min) clusters. If there is a tie, just use the method whose column comes first
'''
def getMinMethod(cluster_table):
    min = len(cluster_table) # start with the min at the number of cells
    min_method = '' # the method with the fewest clusters

    # iterate through methods 
    for method in cluster_table.columns:
        num_clusters = len(pd.unique(cluster_table[method])) # calucate the number of clusters
        if num_clusters < min: # if min, update values
            min_method = method
            min = num_clusters

    return min_method


'''
Construct a confusion matrix from clusters from different methods.
'''
def getMatrix(min_clusters_table, m_clusters_table):
    min_clusters = pd.unique(min_clusters_table) # get a list of all cluster labels from the min method (min)
    m_clusters = pd.unique(m_clusters_table) # get a list of all cluster labels from the other method (m)

    matrix = {} # save the matrix as a dict
    for c1 in min_clusters:
        for c2 in m_clusters:
            set1 = set(min_clusters_table[min_clusters_table == c1].index)
            set2 = set(m_clusters_table[m_clusters_table == c2].index)
            matrix[(c1,c2)] = jaccard_index(set1, set2)

    return matrix


'''
Get clustering method re-labeling key
'''
def getNewLabels(min_clusters_table, m_clusters_table):
    label_key = {} # a dict where the key is the cluster label in the 'm' method and the value is the cluster label in the 'min' method that has the highest jaccard index

    # lists of cluster labels in each method
    m_clusters = pd.unique(m_clusters_table)
    min_clusters = pd.unique(min_clusters_table)

    # for every cluster label in the method, find the most similar cluster label in the min method (using jaccard index)
    for c_label1 in m_clusters:
        indices = {}
        for c_label2 in min_clusters:
            c_label1_items = set(m_clusters_table[m_clusters_table == c_label1].index)
            c_label2_items = set(min_clusters_table[min_clusters_table == c_label2].index)
            i = jaccard_index(c_label1_items, c_label2_items)
            indices[i] = c_label2
        max_index = max(indices.keys()) # find the largest jaccard index
        label_match =  indices[max_index] # get the label in the min method with that max jaccard index
        label_key[c_label1] = label_match # add to key dict
    
    return label_key


'''
Get the cluster label that has the majority vote
'''
def vote(labels):
    c = Counter(labels) # count the most common label

    # check for majority
    label, count = c.most_common()[0] # get most common label and its count
    freq = count / len(labels) # the frequency of the most common label

    if freq > 0.5:
        return label # return mode label if it has majority
    else:
        return np.nan # if it is a tie return NaN


'''
Consensus cluster.
'''
def consensusCluster(cluster_table):

    # compare all methods against the one with the least number of clusters
    min_method = getMinMethod(cluster_table) # get the method with the fewest number of clusters to compare all others against
    methods = list(cluster_table.columns) # get list of methods
    methods.remove(min_method) # remove min method because we will not compare it against itself

    for m in methods:
        label_key = getNewLabels(cluster_table[min_method], cluster_table[m]) # a dict where the key is the cluster label in the 'm' method and the value is the cluster label in the 'min' method that has the highest jaccard index
        # cluster_table[m].replace(label_key, inplace=True) # relabel clusters SLOWER
        cluster_table[m] = cluster_table[m].map(label_key).fillna(cluster_table[m]) # relabel clusters FASTER

    # vote
    consensus_labels = [vote(row) for row in cluster_table.to_numpy()] # a list of clusters labels (or nan) ordered for each item ID in cluster_table
    cluster_table[CONSENSUS] = consensus_labels # add results to table
    return cluster_table

    # make dict with key = index (aka item ID) and value is consensus cluster label
    # consensus_results = dict(zip(cluster_table.index, consensus_labels))
    # return consensus_results
    


'''
Main.
'''
if __name__ == '__main__':

    # constants
    CLUSTER = 'Cluster' # the header of the cluster assignmnet column
    METHOD = 'Method' # the header of the method column
    CONSENSUS = 'Consensus' # the header of the consensus cluster label

    # parse arguments
    args = parseArgs()

    id = args.id # the header of the sample/cell/drug/item ID column

    # get user-defined output dir (strip last slash if present) or set to current
    if args.output is None:
        output = '.'
    elif args.output[-1] == '/':
        output = args.output[:-1]
    else:
        output = args.output


    # if input file has correct format, make a matrix from it
    cluster_table = pd.DataFrame() # a df where cell ID is the index and each column is the cluster assignment from a different algorithm
    for file in args.input:
        if validInput(file): # check that input file has the needed columns
            data, method = getData(file)
            cluster_table[method] = data[CLUSTER]
        else:
            print(f'{file} is incorrectly formatted.')

    # res = consensusCluster(cluster_table)
    # cluster_table[CONSENSUS] = cluster_table.index
    # cluster_table[CONSENSUS] = cluster_table[CONSENSUS].map(res).fillna(cluster_table[CONSENSUS])

    sub_cluster_table = consensusCluster(cluster_table) 
    cluster_table.loc[sub_cluster_table.index] = sub_cluster_table
    # cluster_table.combine_first(sub_cluster_table) # ALTERNATIVE METHOD. NOT SURE WHICH IS FASTER

    print(cluster_table)