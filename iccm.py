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
from sklearn.metrics import jaccard_score

'''
Parse arguments.
'''
def parseArgs():
    parser = argparse.ArgumentParser(description='Implementation of the consensus clustering method proposed in https://www.sciencedirect.com/science/article/pii/S131915781930597X.')
    parser.add_argument('-i', '--input', help='Input cell cluster assignment files.', nargs='*', action="store", dest="input", required=True)
    parser.add_argument('-o', '--output', help='The directory to which output files will be saved.', type=str, required=False)
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
    NEEDED_COLUMNS = [CELLID, CLUSTER, METHOD]

    input_columns = getHeader(file) # get input file header columns
    valid = all(col in input_columns for col in NEEDED_COLUMNS) # check if all needed columns are present in input columns list
    
    return valid


'''
Read file into dataframe.
'''
def getData(file):

    data = pd.read_csv(file, delimiter=',', index_col=CELLID) # load data from input csv into dataframe
    method = data[METHOD].iloc[0] # get method name from first row (this is assuming at all rows are from the same method) NOTE: we may not want to assume this!!
    data = data.drop(METHOD, axis=1) # drop method column
    
    return data, method


'''
Compute all pairwise comparisons of clusters from two methods
'''
def getComparisons(clusters1, clusters2):
    pass


'''
Compute the jaccard index of two clusters (sets)
'''
def jaccard(cluster1, cluster2):
    # return jaccard_score(cluster1, cluster2)
    # OR
    # i = len(cluster1.intersection(cluster2))
    # u = len(cluster1.union(cluster2))
    # jaccard_index = i / u
    # return jaccard_index
    pass





'''
Main.
'''
if __name__ == '__main__':

    # constants
    CELLID = 'CellID' # the header of the cell ID column
    CLUSTER = 'Cluster' # the header of the cluster assignmnet column
    METHOD = 'Method' # the header of the method column

    # parse arguments
    args = parseArgs()

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

    