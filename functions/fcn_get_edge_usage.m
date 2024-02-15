function [edge_usage,percent_usage] = fcn_get_edge_usage(L)
%FCN_GET_EDGE_USAGE     calculate edge usage in shortest paths
%
%   [edge_usage,percent_usage] = fcn_get_edge_usage(L) takes as input a
%   cost matrix and returns the number of times an edge {i,j} was
%   traversed across all shortest paths and the fraction of all edges used
%   in at least one shortest path.
%
%   Inputs:
%               L,     [node x node] cost matrix where weights are already
%                      transformed from measures of affinity to cost.
%   Outputs:
%     
%      edge_usage,    [node x node] matrix of edge counts across all
%                     shortest paths.
%   percent_usage,    [node x node] matrix that transforms usage counts to
%                     percent usage.
%                   
%   NOTE: to calculate shortest paths edge weights must be transformed from
%   measures of affinity to cost. There are many strategies for doing so,
%   including:
%
%       COST = WEIGHT^gamma                 gamma < 0
%       COST = -log(WEIGHT/WEIGHT_MAX)      WEIGHT_MAX = strongest edge
%
%   In Tanner et al, we adopted the first option with gamma = -1.
%
%   NOTE: shortest paths are not well-defined for signed networks. In
%   Tanner et al, we addressed this by transforming the signed matrix so
%   that its entries were all positive. We then adopted the same
%   weight-to-cost transformation as above. Specifically, we performed the
%   following steps:
%
%       1. calcluate minimum nonzero entries in B:
%           Bmin = min(nonzeros(B(:)));
%
%       2. calculate range of nonzero entries in B:
%           Brange = range(nonzeros(B))
%
%       3. subtract the minimum weight from all nonzero entries:
%           B(B ~= 0) = B(B ~= 0) - Bmin.
%
%       4. add a small value to each edge so that the minimum weight is not
%          exactly zero.
%           scale_factor = 0.0001;
%           B(sc ~= 0) = B(sc ~= 0) + Brange*scale_factor;
%
%       5. calculate the cost matrix
%           L = 1./B;
%
%   NOTE: This function uses MATLAB's "graphshortestpath", "graph", and
%   "digraph" functions. In principle, the call to "graphshortestpath" can
%   be replaced with any function that returns the sequence of nodes from a
%   source, s, to a target, t.
%
%   Jacob Tanner, Richard Betzel, 2024

%% get shortest path usage

% number of nodes
n = length(L);

% edge usage matrix
edge_usage = zeros(n);

[~,hops,pmat] = distance_wei_floyd(L);

% loop over source nodes
for s = 1:n
    
    % loop over target nodes
    for t = 1:n
        
        % exclude self
        if s ~= t
            % for source/target pair, return shortest path
            pth = retrieve_shortest_path(s,t,hops,pmat);
            
            % loop over all hops
            for j = 1:length(pth) - 1
                
                % update edge usage matrix
                edge_usage(pth(j),pth(j + 1)) = ...
                    edge_usage(pth(j),pth(j + 1)) + 1;
            end
        end
    end
end

% calculate total number of edges in the matrix
total_edges = nnz(~isnan(L) & ~isinf(L));

% number of edges used at least once
used = nnz(edge_usage);

% percent usage
percent_usage = used/total_edges;