function [B,tspred,tsorig,corr_pred_obs,MSE,local_error] = fcn_fit_model(sc,ts)
%fcn_fit_model    reweight fiber tracts using functional data
%
%   [B,tspred,corr_pred_obs,MSE, ocal_error] = fcn_fit_model(sc,ts)
%   uses a matrix of interregional coupling (sc) to constrain regression
%   model of brain activity. The activity of region i at time (t + 1) is
%   explained by the history of i's connected neighbors at time (t).
%
%   Inputs:
%              sc,    [node x node] coupling matrix
%              ts,    [time x node] matrix of nodal activity for single
%                     scan. can fit using multiple scans by storing each
%                     [time x node] matrix in a cell array whose length is
%                     equal to the number of scans.
%
%   Outputs:
%     
%               B,   [node x node] matrix of regression coefficients; each
%                    non-zero entry in sc is assigned a regression weight.
%          tspred,   regional time series predicted by the model
%   corr_pred_obs,   correlation between the the predicted and observed time
%                    series
%             MSE,   mean squared error between the predicted and observed 
%                    time series
%     local_error,   measures of regional performance (column 1 is MSE and
%                    column 2 is correlation)
%
% Reference:
%   Tanner, J., Faskowitz, J., Teixeira, A. S., Seguin, C., Coletta, L., 
%   Gozzi, A., ... & Betzel, R. F. (2022). Redefining the connectome: A 
%   multi-modal, asymmetric, weighted, and signed description of anatomical 
%   connectivity. bioRxiv, 2022-12.
%
% Jacob Tanner, Richard Betzel, Indiana University, 2024

%% zscore time series
if iscell(ts)
    z = cell(size(ts));
    for scan = 1:length(ts)
        z{scan} = zscore(ts{scan});
    end
else
    z{1} = zscore(ts);
end

%% assign SC 
G = sc;

%% run linear regression

% number nodes
n = length(G);

% to store beta weights
B = zeros(n);

% to store predicted time series
tspred = zeros((size(z{1},1) - 1)*length(z),size(z{1},2));

% loop over nodes
for j = 1:n

    % connected neighbors of node j
    jdx = G(:,j) > 0;

    %jdx(j) = true; % uncommenting this line adds in the diagonal

    % time series at t = 2 to t = T
    y = zeros((size(z{1},1) - 1)*length(z),1);
    x = zeros((size(z{1},1) - 1)*length(z),sum(jdx));

    % loop over scans
    count = 0;
    for scan = 1:length(z)
        yscan = squeeze(z{scan}(2:end,j));
        idx = (count + 1):(count + length(yscan));
        count = count + length(yscan);
        y(idx) = yscan;
        xscan = squeeze(z{scan}(1:end - 1,jdx));
        x(idx,:) = xscan;
    end

    % regressors
    X = [ones(size(x,1),1),x];

    % regression
    b = regress(y,X);

    % store weights
    B(j,jdx) = b(2:end);
    
    % make prediction
    ypred = X*b;

    % store prediction
    tspred(:,j) = ypred;
    
end

%% calculate measures of model performance

% correlation of time series
tsorig = zeros(size(tspred));
count = 0;
for scan = 1:length(z)
    tsscan = z{scan}(2:end,:);
    idx = (count + 1):(count + size(tsscan,1));
    tsorig(idx,:) = tsscan;
    count = count + length(yscan);
end
corr_pred_obs = corr(tspred(:),tsorig(:));

% mean squared error
MSE = nanmean((tspred(:) - tsorig(:)).^2);

% local (node-level) mse
local_error = zeros(n,2);
local_error(:,1) = nanmean((tspred - tsorig).^2);

% local (node-level) correlations
for n = 1:size(tspred,2)
    local_error(n,2) = corr(tspred(:,n),tsorig(:,n));
end



