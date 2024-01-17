function b = fcn_get_geometric_null(sc)


%
%
%   
%
%   Inputs:
%       sc,
%           weighted structural connectome (e.g. AWS connectome, or fiber density connectome)
%       
%              
%
%   Outputs:
%     
%       b,
%           null connectivity matrix
%       
%           
%
%   Example:
%        b = fcn_get_geometric_null(sc);
%        gamma = 1;
%        [M,Q] = community_louvain(b,gamma,[],'negative_asym');

    
    mtrx = sc;
    mask = mtrx ~= 0;
    g = nanmean(mtrx(mask));
    b = mtrx - (mask*g);
