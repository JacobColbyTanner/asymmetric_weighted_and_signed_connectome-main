function b = fcn_get_geographic_null(sc)
%FCN_GET_GEOGRAPHIC_NULL   returns modularity matrix for geographic null
%
%   b = fcn_get_geographic_null(sc) takes an sc matrix and uses the
%   geographic null model (see Bassett et al 2015, Soft Matter or Giusti et
%   al 2016, PRE). This matrix can then serve as input into a modularity
%   maximization function.
%
%   Inputs:
%              sc,     [node x node] connectivity matrix.
%   Outputs:
%
%               b,    [node x node] modularity matrix (sc - expected
%                     weights under geographic null).
%
%   Example:
%        b = fcn_get_geometric_null(sc);
%        gamma = 1;
%        [M,Q] = community_louvain(b,gamma,[],'negative_asym');
%
%   References:
%       Bassett, D. S., Owens, E. T., Porter, M. A., Manning, M. L., & 
%       Daniels, K. E. (2015). Extraction of force-chain network 
%       architecture in granular materials using community detection. Soft 
%       Matter, 11(14), 2731-2744.
%
%       Giusti, C., Papadopoulos, L., Owens, E. T., Daniels, K. E., & 
%       Bassett, D. S. (2016). Topological and geometric measurements of 
%       force-chain structure. Physical Review E, 94(3), 032909.
%
%   Jacob Tanner, Richard Betzel, 2024

% mask of nonzero entries
mask = sc ~= 0;

% mean weight
g = nanmean(sc(mask));

% observed - mean weight
b = sc - (double(mask)*g);
