function [inout_sim,p] = fcn_get_in_out_similarity(B)

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
%       inout_sim,
%           an array of similarity values for each node indicating the similarity of the in_weights and out_weights in the AWS connectome
%       p,
%           p-value for each similarity value
%           
%
%   Example:
%        [B,tspred,corr_pred_obs,MSE, local_error] = fcn_run_sc_regress(sc,ts);
%        [inout_sim,p] = fcn_get_in_out_similarity(B)

    id = B == 0;
    B(id) = nan;
    
    for i = 1:size(B,1)
        
        out = B(i,:)';
        inn = B(:,i);

        [r(i),p(i)] = corr(out,inn,"rows","complete");

    end
    

    inout_sim = r;


