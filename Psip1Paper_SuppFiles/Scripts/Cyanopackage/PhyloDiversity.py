import pickle

import numpy as np
from skbio.stats.distance import mantel


def mantel_pairwise(inputfile):
    global result
    result = list([])
    multi_array_gene = pickle.load(open(inputfile, 'rb'))

    for array in multi_array_gene:

        gene1, gene2, subtree1, subtree2, shared_taxa = array

        if len(shared_taxa) >= 3:
            mantel_test = mantel(subtree1, subtree2, method='pearson', permutations=999)
            if mantel_test[1] < 0.05:
                result.append([[gene1, gene2, 1 - ((mantel_test[0] - (-1)) / (1 - (-1)))],
                               [gene2, gene1, 1 - ((mantel_test[0] - (-1)) / (1 - (-1)))]])

        else:
            if len(shared_taxa) == 0:
                result.append([[gene1, gene2, 1], [gene2, gene1, 1]])

            else:
                result.append([[gene1, gene2, (1 - (
                        (2 * len(shared_taxa)) / (len(subtree1.columns) + len(subtree2.columns))))],
                               [gene2, gene1, (1 - (
                                       (2 * len(shared_taxa)) / (len(subtree1.columns) + len(subtree2.columns))))]])
    return result


class PhyloDiversity:
    """ Aglomerative Hierarchical Clustering based on individual gene trees"""

    def __init__(self, myList: list, df_database: dict, taxa_database: dict, df_storage: np.array):
        self.myList = myList
        self.df_database = df_database
        self.taxa_database = taxa_database
        self.df_storage = df_storage

    """Second part will run a for loop to start producing the pairwise comparison,
     producing the loop to generate data at each pairwise."""

    """Second part will run a for loop to start producing the pairwise comparison,
         producing the loop to generate data at each pairwise."""

    def extraction_matrix(self, reference: str, array_collector: np.array, uniqueList: list):

        gene_data1 = self.df_database[reference]  # Reference dataframe
        taxa_gen1 = self.taxa_database[reference]

        for y in self.myList:

            if y not in uniqueList:
                gene_data2 = self.df_database[y]
                taxa_gen2 = self.taxa_database[y]

                intersec = list(set(taxa_gen1) & set(taxa_gen2))

                pruned_df1 = gene_data1.loc[intersec, intersec]
                pruned_df2 = gene_data2.loc[intersec, intersec]

                test_array = np.array([reference, str(y), pruned_df1, pruned_df2, intersec], dtype=object)
                array_collector = np.append(array_collector, [test_array], axis=0)

        uniqueList.append(reference)

        return array_collector, uniqueList

    """ Final estimation of the mantel distance and storage of the data in the empty data frame 
    which will be return at the end of the function. The input file is the result from the previous script"""
