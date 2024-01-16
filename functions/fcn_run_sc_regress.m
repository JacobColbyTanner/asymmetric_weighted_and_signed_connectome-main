function [B,tspred,corr_pred_obs,MSE, local_error] = fcn_run_sc_regress(sc,ts)

   %Function to get an asymmetric, weighted and signed (AWS) connectome from a binary structural network and activity time series
%
%   This function is a fast and simple way to re-weight structural connectome data using functional imaging data
%
%   Inputs:
%       sc,
%           directed/undirected adjacency matrix describing the anatomical connections between nodes
%       ts,
%           activity time series representing the activity of the same nodes found in the structural connectivity matrix (sc) [time x node]
%              
%
%   Outputs:
%     
%       B,
%           the coefficients from the linear regression model that represent the new weights in the AWS connectome
%            
%       tspred,
%           the predicted time series using the model
%
%       corr_pred_obs,
%           correlation between the the predicted and observed time series
%        
%       MSE, 
%           mean squared error between the predicted and observed time series
%
%       local_error,
%           cell array of the MSE and corr_pred_obs per brain region/node
%
%   Example:
%       [B,tspred,corr_pred_obs,MSE, local_error] = fcn_run_sc_regress(sc,ts)
%       
%
%   References:
%       Tanner et al (2022). Reweighting the connectome



%% zscore time series
tss = zscore(ts);

%% assign SC 
G = sc;
%% run linear regression

n = length(G);

% to store beta weights
B = zeros(n);

% to store predicted time series
tspred = zeros(size(tss));
tspred(1,:,:) = [];

% loop over nodes
for j = 1:n
    jdx = G(:,j) > 0;

    %jdx(j) = true; % uncommenting this line adds in the diagonal

    % the thing we want to predict
    y = squeeze(tss(2:end,j,:));
    y = reshape(y,numel(y),1);

    % the activity of neighbors we'll use for prediction
    x = tss(1:end - 1,jdx,:);
    X = [];
    for i = 1:size(x,3)
        X = [X; x(:,:,i)];
    end
    X = [ones(size(X,1),1),X];

    % regression
    b = regress(y,X);

    % store weights
    B(j,jdx) = b(2:end);



    % make prediction
    ypred = X*b;

    % store prediction
    tspred(:,j,:) = reshape(ypred,[length(ypred)/size(tss,3),size(tss,3)]);

    
    
end

tmp = ts(2:end,:,:);
corr_pred_obs = corr(tspred(:),tmp(:),"rows","complete");

MSE = nanmean((tspred(:)-tmp(:)).^2,"all");


pred_FC = corr(tspred);

act_FC = corr(tmp);


local_error.FC.predicted_FC = pred_FC;
local_error.FC.actual_FC = act_FC;



local_MSE = nanmean((tspred-tmp).^2);

local_error.local_MSE = local_MSE;


for n = 1:size(tspred,2)
    local_corr_pred_obs(n) = corr(tspred(:,n),tmp(:,n),"rows","complete");
    
end

local_error.local_corr_pred_obs = local_corr_pred_obs;



