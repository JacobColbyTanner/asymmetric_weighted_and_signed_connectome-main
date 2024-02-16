# Asymmetric Weighted and Signed Connectome toolbox

This toolbox presents a set of functions that can be used to create and analyze asymmetric, weighted and signed anatomical networks as presented in the paper:
Tanner, J., Faskowitz, J., Teixeira, A. S., Seguin, C., Coletta, L., Gozzi, A., ... & Betzel, R. (2022). Reweighting the connectome: A multi-modal, asymmetric, weighted, and signed description of anatomical connectivity.

Firstly, in order to create such a network you must have two things:

  1. A binarized structural brain network (directed or undirected) describing the anatomical connections between all nodes in the network. This network should be in     the form of an adjacency matrix.

  2. Brain activity time series with the same number of nodes as the binarized structural brain network. This time series data should be organized into a matrix where   the rows represent time points and the columns represent node activity (time x node).

The main.m file uses some example data from a de-identified subject from Human Connectome Project dataset in order to create and plot an asymmetric, weighted and signed network from (1) a binarized structural brain network, and (2) a time series of brain activity. Additionally, this code will also call a variety of functions to analyze this new network and compare it with a fiber density weighted network, the results of which are plotted.

The functions included in this toolbox are listed below. More detailed comments and instructions for usage can be found within each function.

  1. <code>fcn_fit_model</code> : takes in a binary structural network and brain activity time series and returns an asymmetric, weighted and signed connectome, as well as various performance metrics.

  2. <code>fcn_get_asymmetry</code> : returns three measures of asymmetry. The asymmetry matrix where each entry is the difference between the weight of edge(i,j) and edge(j,i), the absolute asymmetry matrix where each entry is the absolute difference between edge(i,j) and edge(j,i), and the sign asymmetry matrix that describes where edge(i,j) and edge(j,i) had different signs.

  3. <code>fcn_get_in_out_similarity</code> : returns an array of values (one for each node in the network) that describe the similarity of the in-weights and out-weights for each node. The in-weights are the weights that were used by a linear regression model to predict the future activity of this node. So, the in-weights are weights describing when other nodes are used to predict this nodes activity. In contrast, the out-weights are when this nodes activity is used to predict a different nodes activity.

  4. <code>fcn_consensus_communities</code> : takes in multiple partitions of the same connectivity matrix into communities/modules and returns a consensus partition. This function includes the use of two functions to relabel and identify unique partitions (<code>fcn_relabel_partitions</code> and <code>fcn_unique_partitions</code> respectively).

  5. <code>fcn_sort_communities</code> : relabels the communities of a partition such that the largest community gets labeled 1, the next largest gets labeled 2, and so on.

  6. <code>fcn_get_geometric_null</code> : creates a null connectivity matrix for usage in modularity maximization. This null assigns every existing edge in the network to the mean weight value across all existing edges.

  7. <code>fcn_bilaterality</code> : takes in a set of community labels as well as labels for the hemisphere each node belongs to and computes a measure of laterality. This measure describes how much each of the modules tends to be concentrated in one of the two hemispheres.

  8. <code>fcn_get_edge_usage</code> : counts the number of times each edge is used in the shortest paths between all nodes in the network.

In addition, this toolbox uses a number of functions from the Brain Connectivity toolbox which is publicly available at:
https://sites.google.com/site/bctnet/.

If you use this toolbox in your published work, please cite our paper:

  - Tanner, Jacob, Joshua Faskowitz, Andreia Sofia Teixeira, Caio Seguin, Ludovico Coletta, Alessandro Gozzi, Bratislav Mišić, and Richard F. Betzel. "Redefining the      connectome: A multi-modal, asymmetric, weighted, and signed description of anatomical connectivity." bioRxiv (2022): 2022-12.    
    https://doi.org/10.1101/2022.12.19.519033.








