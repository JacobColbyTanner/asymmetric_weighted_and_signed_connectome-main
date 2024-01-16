# asymmetric_weighted_and_signed_connectome

This toolbox presents a set of functions that can be used to create and analyze asymmetric, weighted and signed anatomical networks as presented in the paper:

Tanner, J., Faskowitz, J., Teixeira, A. S., Seguin, C., Coletta, L., Gozzi, A., ... & Betzel, R. (2022). Reweighting the connectome: A multi-modal, asymmetric, weighted, and signed description of anatomical connectivity. 

Firstly, in order to create such a network you must have two things:

(1) A binarized structural brain network (directed or undirected) describing the anatomical connections between all nodes in the network. This network should be in the form of an adjacency matrix.

(2) Brain activity time series with the same number of nodes as the binarized structural brain network. This time series data should be organized into a matrix where the rows represent time points and the columns represent node activity (time x node).

The main_test.m file uses some example data from a de-identified subject from Human Connectome Project data in order to create and plot an asymmetric, weighted and signed network from (1) a binarized structural brain network, and (2) a time series of brain activity. Additionally, this code will also call a variety of functions to analyze this new network and compare it with a fiber density weighted network, the results of which are plotted. 

The analysis functions included in this toolbox consist of:

Asymmetry analysis functions: these functions analyze the asymmetries in these new networks by comparing the weights of edges from node i to node j, with the edges from node j to node i.

  <code> fcn_get_asymmetry <code/>: returns a matrix where each entry is the difference between the weight of edge(i,j) and edge(j,i)

  fcn_get_sign_asymmetry: returns a binary matrix where entries of value 1 indicate that edges (i,j) and (j,i) had a different sign(+/-).

  fcn_get_in_out_similarity: returns an array of values (one for each node in the network) that describe the similarity of the in-weights and out-weights for each node. The in-weights are the weights that were used by a linear regression model to predict the future activity of this node. So, the in-weights are weights describing when other nodes are used to predict this nodes activity. In contrast, the out-weights are when this nodes activity is used to predict a different nodes activity.

  community_louvain: this is a function from the brain connectivity toolbox (https://sites.google.com/site/bctnet/) that computes modularity maximization with the Louvain algorithm in order to divide an adjacency matrix into modules. For more info on modularity maximization, see: Esfahlani, F. Z., Jo, Y., Puxeddu, M. G., Merritt, H., Tanner, J. C., Greenwell, S., ... & Betzel, R. F. (2021). Modularity maximization as a flexible and generic framework for brain network exploratory analysis. Neuroimage, 244, 118607.

  get_geometric_null: creates a null connectivity matrix for usage in modularity maximization. This null assigns every existing edge in the network to the mean weight value across all existing edges.

  fcn_bilaterality: takes in a set of community labels (produced by community_louvain) as well as labels for the hemisphere each node belongs to and computes a measure of laterality. This measure describes how much each of the modules tends to be concentrated in one of the two hemispheres.

  get_edge_usage: counts the number of times each edge is used in the shortest paths between all nodes in the network.

  randmio_und: this is a function from the brain connectivity toolbox (https://sites.google.com/site/bctnet/) that randomizes an undirected network while preserving the degree distribution.


If you use this toolbox in your published work, please cite our paper:

Tanner, J., Faskowitz, J., Teixeira, A. S., Seguin, C., Coletta, L., Gozzi, A., ... & Betzel, R. (2022). Reweighting the connectome: A multi-modal, asymmetric, weighted, and signed description of anatomical connectivity.









