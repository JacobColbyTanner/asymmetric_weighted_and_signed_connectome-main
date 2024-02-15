function [r] = fcn_get_in_out_similarity(mat)
%FCN_GET_IN_OUT_SIMILARITY      similarity of input/output weights
%
%   [r] = fcn_get_in_out_similarity(mat) calculates the similarity of 
%   input/output weights for a connectivity matrix, mat, as a correlation 
%   coefficient.
%
%   Inputs:
%             mat,     [node x node] connectivity matrix.
%
%   Outputs:
%     
%               r,     [1 x node] vector of similarity
%
%   Example:
%
%   >> B = fcn_run_sc_regress(sc,ts);
%   >> r = fcn_get_in_out_similarity(B);
%
%   Jacob Tanner, Richard Betzel, 2024

%% calculate in-out similarity
% number of nodes
n = length(mat);

% zscore columns
zc = zscore(mat);

% zscore rows
zr = zscore(mat');

% calculate similarity as correlation coefficient
r = sum(zc.*zr)/(n - 1);