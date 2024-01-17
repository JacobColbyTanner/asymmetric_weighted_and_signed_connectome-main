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
%        b = fcn_get_geometric_null;
%        [edge_usage,percent_usage] = fcn_get_edge_usage(B)

    
    mtrx = sc;
    mask = mtrx ~= 0;
    g = nanmean(mtrx(mask));
    b = mtrx - (mask*g);
