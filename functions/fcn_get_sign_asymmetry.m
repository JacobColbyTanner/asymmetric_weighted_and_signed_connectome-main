function [sign_asym] = fcn_get_sign_asymmetry(B)


%
%
%   
%
%   Inputs:
%       B,
%           the coefficients from the linear regression model that represent the new weights in the AWS connectome
%       
%              
%
%   Outputs:
%     
%       sign_asym,
%           a binary matrix where each 1 indicates where edge(i,j) and edge(j,i) have a different sign (+/-) in the AWS connectome
%
%           
%
%   Example:
%        [B,tspred,corr_pred_obs,MSE, local_error] = fcn_run_sc_regress(sc,ts);
%        [sign_asym] = fcn_get_sign_asymmetry(B)



    pos = B > 0;
    neg = B < 0;
    %get elements for i,j is pos/neg and j,i is neg,pos
    sign_asym = pos+neg' == 2;
