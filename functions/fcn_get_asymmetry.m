function [asymmetry, abs_asymmetry] = fcn_get_asymmetry(B)

%Function to measure the laterality of a modular partition. Laterality measures how much a module favors one hemisphere over the other.
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
%       asymmetry,
%           asymmetry matrix where each value gives the difference between edge(i,j) and edge(j,i)
%       abs_asymmetry,
%           absolute asymmetry matrix where each value gives the absolute difference between edge(i,j) and edge(j,i)
%           
%
%   Example:
%        [B,tspred,corr_pred_obs,MSE, local_error] = fcn_run_sc_regress(sc,ts);
%        [asymmetry, abs_asymmetry] = fcn_get_asymmetry(B);





    asymmetry = B-B';
    abs_asymmetry = abs(B-B');

