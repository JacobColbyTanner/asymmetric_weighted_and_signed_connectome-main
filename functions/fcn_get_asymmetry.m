function [asymmetry,abs_asymmetry,sign_sym] = fcn_get_asymmetry(B)
%FCN_GET_ASYMMETRY    returns weight/sign asymmetries
%
%   [asymmetry,abs_asymmetry,sign_sym] = fcn_get_asymmetry(B) takes as 
%   input the matrix of regression coefficients and returns, as output, two
%   matrices whose elements correspond to differences in in-weight versus
%   out-weight (and the absolute value).
%   
%   Inputs:
%             B,     [node x node] matrix of regression coefficients; each
%                    non-zero entry in sc is assigned a regression weight.
%
%   Outputs:
%     
%      asymmetry,    [node x node] matrix whose element {i,j} is equal to
%                    B(i,j) - B(j,i).
%  abs_asymmetry,    same as asymmetry matrix, but with elements equal to
%                    |B(i,j) - B(j,i)| where |x| denotes the absolute value
%                    of x.
%       sign_sym,    [node x node x 3] matrix that returns asymmetries in
%                    the sign of connection weights:
%                       1. e_{ij} > 0 and e_{ji} > 0
%                       2. e_{ij} < 0 and e_{ji} < 0
%                       3. sign(e_{ij}) ~= sign(e_{ji})
%
% Jacob Tanner, Richard Betzel, Indiana University, 2024
%

%% calculate differences between matrix and its transpose
asymmetry = B - B';
abs_asymmetry = abs(asymmetry);

Bsign = sign(B);
sign_sym = zeros([size(B),3]);
sign_sym(:,:,1) = (Bsign > 0) & (Bsign' > 0);
sign_sym(:,:,2) = (Bsign < 0) & (Bsign' < 0);
sign_sym(:,:,3) = (Bsign > 0) & (Bsign' < 0);
sign_sym(:,:,3) = sign_sym(:,:,3) | sign_sym(:,:,3)';