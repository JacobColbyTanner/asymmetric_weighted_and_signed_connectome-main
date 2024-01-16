
%%Load data and add path to functions___________________________________
addpath("functions/")
load("data/sample_data.mat")
load("data/hcp400.mat")
load("data/hcp400centroids.mat")
load("data/schaefer-yeo17_400node_permuted_inds.mat")

%%Run AWS model and plot weight matrices___________________________________
[B,tspred,corr_pred_obs,MSE, local_error] = fcn_run_sc_regress(sc,ts);
save("data/AWS_connectome.mat","B")


figure;
subplot(1,2,1)
imagesc(sc)
title("Original fiber-density weighted connectome")
axis("square")
colorbar
subplot(1,2,2)
imagesc(B,[-1 1])
title("Asymmetric, weighted, and signed connectome")
axis("square")
colorbar


%%calculate asymmetry values___________________________________
%regular and absolute asymmetry
[asymmetry, abs_asymmetry] = fcn_get_asymmetry(B);
figure;
subplot(1,2,1)
imagesc(asymmetry, [-1 1])
title("Asymmetry matrix")
axis("square")
colorbar
subplot(1,2,2)
imagesc(abs_asymmetry)
title("Absolute asymmetry matrix")
axis("square")
colorbar


%sign asymmetry: calculate where edge(i,j) and edge(j,i) have a different
%sign
[sign_asym] = fcn_get_sign_asymmetry(B);
figure;
imagesc(sign_asym,[0 1])
title("Sign asymmetric edges")
axis("square")
colorbar

%in_out similarity: calculate the correlation between in_weights and
%out_weights per region
load("data/hcp400.mat")

[inout_sim,p] = fcn_get_in_out_similarity(B);

systems = nan(size(B,1),max(lab));
for i = 1:max(lab)
    idx = lab == i;
    systems(1:sum(idx),i) = inout_sim(idx);
    
end

figure;
boxplot(systems)
xticks(1:16)
xticklabels(net)
ylabel("in/out similarity")







%%calculate modules and laterality___________________________________

%choose null model for community detection (uncomment to select)
null = "geo" %or null = "neg_asym"

if null == "geo"
    null_model =  get_geometric_null(B); %custom geometric null model from our paper
else
    null_model = 'negative_asym'; %neg asymmetric null model
end
num_iter = 10; %change to iteration preference (10 for example purposes only)
lat_AWS = [];
for i = 1:num_iter
    %get modules AWS
    gamma = 1;
    [M,Q] = community_louvain(B,gamma,[],null_model);
    hemi = ones(400,1);
    hemi(201:400) = 2;
    nrand = 100;
    %get laterality AWS
    lat_AWS(i) = fcn_bilaterality(M,hemi,nrand);
end

if null == "geo"
    null_model =  get_geometric_null(sc); %custom geometric null model from our paper
else
    null_model = 'negative_asym'; %neg asymmetric null model
end

lat_sc = [];
for i = 1:num_iter
    %get modules sc
    gamma = 1;
    [M,Q] = community_louvain(sc,gamma,[],null_model);
    hemi = ones(400,1);
    hemi(201:400) = 2;
    nrand = 100;
    %get laterality sc
    lat_sc(i) = fcn_bilaterality(M,hemi,nrand);
end

both = [lat_AWS(:), lat_sc(:)];

figure;
boxplot(both)
xticklabels({"AWS", "fiber density"})
ylabel("laterality")



%% get edge usage on shortest path backbone_______________________________

[edge_usage_fiber,percent_usage_fiber] = get_edge_usage(sc);

[edge_usage_AWS,percent_usage_AWS] = get_edge_usage(B);

figure;
subplot(1,2,1)
imagesc(edge_usage_AWS, [0 5])
axis("square")
title("edge usage AWS")
colorbar

subplot(1,2,2)
imagesc(edge_usage_fiber, [0 5])
axis("square")
title("edge usage fiber density")
colorbar

both = [percent_usage_AWS,percent_usage_fiber];

figure;
bar(both)
xticklabels({"AWS", "fiber density"})
ylabel("percent of edges used")
title("shortest paths backbone")

%%Get connection null models & performance________________________________

n = size(sc,1);

%minimally wired null model in which only the shortest (least-costly) connections are preserved 
% (while preserving an equal number of connections)
dist = squareform(pdist(coor));
neworder = randperm(n);
D = dist;
m = nnz(sc)/2;
edges = find(triu(ones(n),1));
Dvec = D(edges);
[~,idxsort] = sort(Dvec,'ascend');
SCneworder = zeros(n);
SCneworder(edges(idxsort(1:m))) = 1;
minwire_sc = SCneworder + SCneworder';


% a reordered null model in which nodes orders were randomly permuted, 
rp = randperm(n);
reorder_sc = sc(rp,:);

% a "spin” model in which nodes orders were randomly permuted but geometry preserved
ri = randi(50000);
ind = PERMS(:,ri);
spin_sc = sc(ind,:);

% a topological null model in which nodes  degrees
numiter = 32;
topological_sc = randmio_und(sc,numiter);

%calculate and plot performances
[BB,tspred,corr_pred_obs,MSE, local_error] = fcn_run_sc_regress(sc,ts);
[BB,tspred,corr_pred_obs_minwire,MSE_minwire, local_error] = fcn_run_sc_regress(minwire_sc,ts);
[BB,tspred,corr_pred_obs_reorder,MSE_reorder, local_error] = fcn_run_sc_regress(reorder_sc,ts);
[BB,tspred,corr_pred_obs_spin,MSE_spin, local_error] = fcn_run_sc_regress(spin_sc,ts);
[BB,tspred,corr_pred_obs_topo,MSE_topo, local_error] = fcn_run_sc_regress(topological_sc,ts);

all_corr = [corr_pred_obs,corr_pred_obs_minwire,corr_pred_obs_reorder,corr_pred_obs_spin,corr_pred_obs_topo];
all_MSE = [MSE,MSE_minwire,MSE_reorder,MSE_spin, MSE_topo];
figure;
subplot(1,2,1)
bar(all_corr)
ylabel("correlation pred vs. obs")
xticks(1:5)
xticklabels({"observed","min rewire","reorder","spin","degree preserving"})

subplot(1,2,2)
bar(all_MSE)
ylabel("MSE")
xticks(1:5)
xticklabels({"observed","min rewire","reorder","spin","degree preserving"})



