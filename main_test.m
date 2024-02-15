clear all
close all
clc

%% Before you run any code

% The script implements many of the analyses from Tanner et al. That
% manuscript contributes a simple method for estimating weights on
% already-defined edges in a structural connectome. Specifically, we use
% linear regression to explain the future activity of node i based on the
% history of its connected neighbors. The regression coefficients obtained
% from the model are treated as the weights from each of i's neighbors to
% i.

% The resulting matrix, which we denote as "B" in the following script, is
% signed and asymmetric. The following sections explore some of its
% properties that were covered in the manuscript.
%
%   1. calculate and plot asymmetries in the connection weights as well as
%      consensus/dissensus in their sign (fcn_get_asymmetry)
%   2. calculate and plot similarity of incoming/outgoing connectivity
%      profiles for each brain region (fcn_get_in_out_similarity.m)
%   3. calculate the community structure of the B matrix and compare its
%      laterality with that of the fiber density matrix (undirected SC)
%      (fcn_bilaterality; fcn_get_geographic_null)
%   4. calculate shortest paths in the signed matrix and compare edge usage
%      with that of the fiber density matrix (fcn_get_edge_usage)
%   5. compare model fitness against null models.

%% Load data and add path to functions
addpath(genpath('functions/'));
load('data/participant_data.mat');
load('data/hcp400.mat');
load('data/hcp400centroids.mat');
load('data/hcp400cmapOther.mat');
load('data/schaefer-yeo17_400node_permuted_inds.mat');

%% Run model and plot matrices

% fit edge weights to coupling mask
[B,tspred,tsorig,corr_pred_obs,MSE,local_error] = fcn_fit_model(sc,ts);

% print fitness
fprintf('\nMean squared error between observed and predicted time series: %.2f\n',MSE);
fprintf('Correlation between observed and predicted time series: %.2f\n',corr_pred_obs);

% plot outcome
figure('position',[1000,1000,1200,120]);

% plot example time series
subplot(1,5,[1,2]);

% select random node and plot empirical and predicted time series
random_node = randi(length(sc));
idxt = 1:size(ts{1}) - 1;
ph = plot(...
    idxt,tspred(idxt,random_node),...
    idxt,tsorig(idxt,random_node));
title(sprintf('region %i',random_node));
legend(ph,{'predicted','original'},'location','northeast')

% plot binarized sc matrix
subplot(1,5,3);
imagesc(sc ~= 0);
title('Coupling mask');
colorbar;

% plot log fiber density matrix
subplot(1,5,4)
imagesc(log10(sc),[-3,0])
title('Fiber density')
colorbar

% weight matrix estimated from model
subplot(1,5,5)
imagesc(B,[-max(abs(B(:))),max(abs(B(:)))]*0.125)
title('Model weights')
colorbar

%% Plot local error statistics

% draw regional mse (how well model fits region i's time series)
% NOTE: smaller values = better fit
figure('position',[1000,1000,900,120]);
subplot(1,3,1);
scatter3(coor(:,1),coor(:,2),coor(:,3),25,local_error(:,1),'filled');
view([-90,90]);
axis image;
colorbar;
title('local MSE');

% draw regional correlations (how well model fits region i's time series)
% NOTE: larger values = better fit
subplot(1,3,2);
scatter3(coor(:,1),coor(:,2),coor(:,3),25,local_error(:,2),'filled');
view([-90,90]);
axis image;
colorbar;
title('local correlations');

% compare
subplot(1,3,3);
scatter(local_error(:,1),local_error(:,2));
xlabel('local MSE'); ylabel('local correlation');
text(max(local_error(:,1)),max(local_error(:,2)),...
    sprintf('spearman rho=%.2f',corr(local_error(:,1),local_error(:,2),'type','spearman')));

%% Calculate edge weight asymmetries

% calculate differences between incoming/outgoing weights for each edge
[asymmetry,abs_asymmetry,sign_asymmetry] = fcn_get_asymmetry(B);

% difference between in/out weight
figure('position',[1000,1000,900,200]);
subplot(1,3,1);
imagesc(asymmetry,[-1,1]*0.25);
title('W_{ij} - W_{ji}')
colorbar

% absolute difference between in/out weight
subplot(1,3,2);
imagesc(abs_asymmetry);
title('|W_{ij} - W_{ji}|')
colorbar

% agreement and disagreement of edge weight sign
mat = sum(bsxfun(@times,sign_asymmetry,permute(3:-1:1,[3,1,2])),3);
subplot(1,3,3);
imagesc(mat);
title('sign asymmetries, 3=++, 2=--, 1=+-')
colorbar;

%% Calculate similarity of incoming and outgoing connectivity profiles

% calculate similarity
r = fcn_get_in_out_similarity(B);

% make boxplot where similarity is grouped by Schaefer systems
figure('position',[1000,1000,700,200]);
subplot(1,2,1);
boxplot(r,lab,'labels',net,'color',cmap)
title('similarity incoming/outgoing connectivity profiles');

% plot same data in anatomical space
subplot(1,2,2);
scatter3(coor(:,1),coor(:,2),coor(:,3),25,r,'filled');
view([-90,90]);
axis image;
colorbar;
title('similarity incoming/outgoing connectivity profiles');

%% Estimate communities and partition laterality

% NOTE: for illustrative purposes we run the community detection a small
% number of times. in principle, the number of iterations should be large.

% number of iterations
numiter = 25;

% preallocate arraays
n = length(B);
ci_sc = zeros(n,numiter);  ci_mdl = ci_sc;             ci_geo = ci_sc;
lat_sc = zeros(numiter,1); lat_mdl = zeros(numiter,1); lat_geo = zeros(numiter,1);

% modularity matrix for geographic model
b = fcn_get_geographic_null(B);

% loop over iterations
for iter = 1:numiter

    % from fiber density matrix
    ci_sc(:,iter) = fcn_sort_communities(community_louvain(sc));

    % from model weights
    ci_mdl(:,iter) = fcn_sort_communities(community_louvain(B,[],[],'negative_asym'));

    % using geographic null model (Fig. S10, S11, and S12 in paper)
    ci_geo(:,iter) = fcn_sort_communities(community_louvain(ones(size(B)),[],[],b));

    % define parameters for calculating laterality
    nrand = 1000;
    hemi = [ones(n/2,1)*1; ones(n/2,1)*2];

    % calculate laterality of each partition
    lat_sc(iter) = fcn_bilaterality(ci_sc(:,iter),hemi,nrand);
    lat_mdl(iter) = fcn_bilaterality(ci_mdl(:,iter),hemi,nrand);
    lat_geo(iter) = fcn_bilaterality(ci_geo(:,iter),hemi,nrand);

end

% plot laterality
figure('position',[1000,1000,200,250]);
boxplot([lat_sc,lat_mdl,lat_geo],'labels',{'fiber density','model weights','geographic'})
ylabel('community laterality, \Lambda');

% calculate consensus partitions
cicon_sc = fcn_consensus_communities(ci_sc,numiter,false);
cicon_mdl = fcn_consensus_communities(ci_mdl,numiter,false);
cicon_geo = fcn_consensus_communities(ci_geo,numiter,false);

% plot communities in anatomical space
figure('position',[1000,1000,900,120]);
subplot(1,3,1);
scatter3(coor(:,1),coor(:,2),coor(:,3),25,cicon_sc,'filled');
view([-90,90]);
axis image;
colorbar;
title('Communities (Fiber Density)');

subplot(1,3,2);
scatter3(coor(:,1),coor(:,2),coor(:,3),25,cicon_mdl,'filled');
view([-90,90]);
axis image;
colorbar;
title('Communities (Model weights)');

subplot(1,3,3);
scatter3(coor(:,1),coor(:,2),coor(:,3),25,cicon_geo,'filled');
view([-90,90]);
axis image;
colorbar;
title('Communities (Geographic null model)');

%% Calculate edge usage and shortest paths backbone

% To this, we must first transform weights from measures of affinity to a
% measure of cost. We adopt the "reciprocal transform".

% weak fibers are now more costly
cost_sc = 1./sc;

% for signed matrix we need to handle this differently
Bcopy = B;
bmin = min(nonzeros(B));
brange = range(nonzeros(B));
Bcopy(B ~= 0) = Bcopy(B ~= 0) - bmin;
scale_factor = 0.0001;
Bcopy(B ~= 0) = Bcopy(B ~= 0) + brange*scale_factor;
cost_b = 1./Bcopy;

% calculate edge usage and shortest paths
[edge_usage_sc,percent_usage_sc] = fcn_get_edge_usage(cost_sc);

% calculate edge usage and shortest paths
[edge_usage_b,percent_usage_b] = fcn_get_edge_usage(cost_b);

% visualize edge usage matrices
cap = 100;
figure('position',[1000,1000,700,200])
subplot(1,2,1)
imagesc(edge_usage_b,[0,cap])
title('Model weights - edge usage')
colorbar;

subplot(1,2,2)
imagesc(edge_usage_sc,[0,cap])
title('Fiber Density - edge usage')
colorbar;

% visualize edge usage
figure('position',[1000,1000,200,250]);
bar(1:2,[percent_usage_sc,percent_usage_b]);
set(gca,...
    'xtick',1:2,...
    'xticklabel',{'fiber density','model weights'})
ylabel({'percentage of edges involved','in at least one shortest path'})

%% Compare model fitness against null models

% In the main text we explored five null models. Here, we provide code for
% generating those models.

% degree-preserving but otherwise random
num_iter_per_edge = 2^5;
sc_degpres = randmio_und(sc,num_iter_per_edge);

% minimally-wired
sc_minwire = zeros(n);
mask = triu(ones(n),1) > 0;
dist = squareform(pdist(coor));
dvec = dist(mask);
[~,idxsort] = sort(dvec,'ascend');
edges = find(mask);
sc_minwire(edges(idxsort(1:(nnz(sc)/2)))) = 1;
sc_minwire = sc_minwire + sc_minwire';

% randomly reordered
r = randperm(n);
sc_randorder = sc(r,r);

% randomly reordered (but spatially constrained)
r = PERMS(:,1);
sc_spin = sc(r,r);

% concatenate matrices
mats = cat(3,sc_degpres,sc_minwire,sc_randorder,sc_spin);

% preallocate arrays for fitness measures
r = zeros(size(mats,3),1); m = r;

% fit models using null sc
for i = 1:size(mats,3)
    [~,~,~,r(i),m(i),~] = fcn_fit_model(mats(:,:,i),ts);
end

% make plots
names = {'degree preserving','min. wire','rand','spin'};

% plot comparison using mse as a measure of fitness
[~,idxsort] = sort(m,'ascend');
figure('position',[1000,1000,200,250]);
bar(m(idxsort)); hold on;
ph = plot([0.5,length(m) + 0.5],MSE*ones(1,2),'r');
legend(ph,'intact data','location','southwest');
set(gca,...
    'xtick',1:length(m),...
    'xticklabel',names(idxsort));
ylabel('mean squared error')

% plot comparison using correlation as a measure of fitness
[~,idxsort] = sort(r,'descend');
figure('position',[1000,1000,200,250]);
bar(r(idxsort)); hold on;
ph = plot([0.5,length(r) + 0.5],corr_pred_obs*ones(1,2),'r');
legend(ph,'intact data','location','southwest');
set(gca,...
    'xtick',1:length(r),...
    'xticklabel',names(idxsort));
ylabel('correlation')