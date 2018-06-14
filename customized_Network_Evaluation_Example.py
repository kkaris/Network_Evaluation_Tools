# This is a customized version of the example here:
# https://github.com/kkaris/Network_Evaluation_Tools/blob/master/Network%20Evaluation%20Examples/Network%20Evaluation%20Example.ipynb

try:
    from network_evaluation_tools import data_import_tools as dit
except ImportError:  # Importing two times solves the ImportError apparently
    from network_evaluation_tools import data_import_tools as dit

from network_evaluation_tools import network_evaluation_functions as nef
from network_evaluation_tools import network_propagation as prop
import pandas as pd
import numpy as np

input_network_file = './interaction_datasets/head3000_stripped_HUMAN_INTERACTOME.sif'
gene_set_file = './Data/DisGeNET_genesets.txt'

# Load network (We choose a smaller network here for the example's sake)
network = dit.load_network_file(input_network_file, verbose=True)

# Load gene sets for analysis
genesets = dit.load_node_sets(gene_set_file)

# Calculate geneset sub-sample rate
genesets_p = nef.calculate_p(network, genesets)

# Determine optimal alpha for network (can also be done automatically by next step)
alpha = prop.calculate_alpha(network)
# print alpha

# Calculate network kernel for propagation
kernel = nef.construct_prop_kernel(network, alpha=alpha, verbose=True)

# Might want to tweak values here to speed up calculation
# Calculate the AUPRC values for each gene set
AUPRC_values = nef.small_network_AUPRC_wrapper(kernel, genesets, genesets_p, n=30, cores=4, verbose=True)

# Construct null networks and calculate the AUPRC of the gene sets of the null networks
# We can use the AUPRC wrapper function for this
null_AUPRCs = []
for i in range(10):
    shuffNet = nef.shuffle_network(network, max_tries_n=10, verbose=True)
    shuffNet_kernel = nef.construct_prop_kernel(shuffNet, alpha=alpha, verbose=False)
    shuffNet_AUPRCs = nef.small_network_AUPRC_wrapper(shuffNet_kernel, genesets, genesets_p, n=30, cores=4, verbose=False)
    null_AUPRCs.append(shuffNet_AUPRCs)
    print 'shuffNet', repr(i+1), 'AUPRCs calculated'

# Construct table of null AUPRCs
null_AUPRCs_table = pd.concat(null_AUPRCs, axis=1)
null_AUPRCs_table.columns = ['shuffNet'+repr(i+1) for i in range(len(null_AUPRCs))]

# Calculate performance metric of gene sets
network_performance = nef.calculate_network_performance_score(AUPRC_values, null_AUPRCs_table, verbose=True)
network_performance.name = 'Test Network'
network_performance.to_csv('./Data/Network_Performance.csv', sep='\t')

# Calculate network performance gain over median null AUPRC
network_perf_gain = nef.calculate_network_performance_gain(AUPRC_values, null_AUPRCs_table, verbose=True)
network_perf_gain.name = 'Test Network'
network_perf_gain.to_csv('./Data/Network_Performance_Gain.csv', sep='\t')

# Rank network on average performance across gene sets vs performance on same gene sets in previous network set
all_network_performance = pd.read_csv('./Data/Network_Performance.csv', index_col=0)
all_network_performance_filt = pd.concat([network_performance, all_network_performance.ix[network_performance.index]], axis=1)
network_performance_rank_table = all_network_performance_filt.rank(axis=1, ascending=False)
network_performance_rankings = network_performance_rank_table['Test Network']

# Rank network on average performance gain across gene sets vs performance gain on same gene sets in previous network set
all_network_perf_gain = pd.read_csv('./Data/Network_Performance_Gain.csv', index_col=0)
all_network_perf_gain_filt = pd.concat([network_perf_gain, all_network_perf_gain.ix[network_perf_gain.index]], axis=1)
network_perf_gain_rank_table = all_network_performance_filt.rank(axis=1, ascending=False)
network_perf_gain_rankings = network_perf_gain_rank_table['Test Network']

# Network Performance
network_performance_metric_ranks = pd.concat([network_performance, network_performance_rankings, network_perf_gain, network_perf_gain_rankings], axis=1)
network_performance_metric_ranks.columns = ['Network Performance', 'Network Performance Rank', 'Network Performance Gain', 'Network Performance Gain Rank']
network_performance_metric_ranks.sort_values(by=['Network Performance Rank', 'Network Performance', 'Network Performance Gain Rank', 'Network Performance Gain'],
                                             ascending=[True, False, True, False])

# Construct network summary table
network_summary = {}
network_summary['Nodes'] = int(len(network.nodes()))
network_summary['Edges'] = int(len(network.edges()))
network_summary['Avg Node Degree'] = np.mean(dict(network.degree()).values())
network_summary['Edge Density'] = 2*network_summary['Edges'] / float((network_summary['Nodes']*(network_summary['Nodes']-1)))
network_summary['Avg Network Performance Rank'] = network_performance_rankings.mean()
network_summary['Avg Network Performance Rank, Rank'] = int(network_performance_rank_table.mean().rank().ix['Test Network'])
network_summary['Avg Network Performance Gain Rank'] = network_perf_gain_rankings.mean()
network_summary['Avg Network Performance Gain Rank, Rank'] = int(network_perf_gain_rank_table.mean().rank().ix['Test Network'])
for item in ['Nodes', 'Edges' ,'Avg Node Degree', 'Edge Density', 'Avg Network Performance Rank', 'Avg Network Performance Rank, Rank',
             'Avg Network Performance Gain Rank', 'Avg Network Performance Gain Rank, Rank']:
    print item+':\t'+repr(network_summary[item])

