function ciu = fcn_consensus_communities(ci,niter,vis)
%FCN_CONSENSUS_COMMUNITIES   generate consensus communities
%
%   ciu = fcn_consensus_communities(ci,niter,vis) takes as input a [node x
%   partition] ensemble, ci, and returns an estimate of consensus
%   communities. The algorithm calculates the coassignment matrix,
%   estimates the expected coassignment, and defines a consensus modularity
%   matrix (coassignment minus expected coassignment). This matrix is
%   clustered using modularity maximization and the entire procedure
%   repeated until convergence.
%
%   Inputs:
%
%       ci,     [node x iteration] matrix of partitions
%    niter,     number of times to cluster co-assignment matrix
%      vis,     true or false (default) for visualizing progress
%      opt,     determines which community detection code to use ('bct' or
%               'genlouvain').
%
%   Outpts:
%
%      ciu,     consensus partitions
%
%   Rick Betzel, Indiana University, 2014
%

if ~exist('vis','var')
    vis = false;
end
n = size(ci,1);
mask = triu(ones(n),1) > 0;
if size(ci,1) == size(ci,2)
    d = ci;
else
    d = fcn_agreement(ci);
end
goFlag = length(unique(d));
if goFlag <= 2
    CiCon = ci;
end

totalreps = 0;
while goFlag > 2

    totalreps = totalreps + 1;

    mu = mean(d(mask));
    b = d - mu;
    b(1:(n + 1):end) = 0;
    CiCon = zeros(n,niter);
    for iRep = 1:niter
        CiCon(:,iRep) = community_louvain(ones(n),[],[],b);
        if vis & mod(iRep,10) == 0
            imagesc(CiCon); drawnow;
        end
    end
    d = agreement(CiCon);
    goFlag = length(unique(d));

    if totalreps >= 10
        [CiUnq,freq] = fcn_unique_partitions(CiCon);
        [~,idxmax] = max(freq);
        CiCon = CiUnq(:,idxmax);
        goFlag = 1;
    end
    
end
ciu = CiCon(:,1);


function [ciu,valu,mult,indx] = fcn_unique_partitions(ci,val)
ci = fcn_relabel_partitions(ci);
if nargin == 1
    val = 1:size(ci,2);
end

ciu = []; mult = []; valu = [];
count = 0;
c = 1:size(ci,2);
while ~isempty(ci)
    count = count + 1;
    tgt = ci(:,1);
    ciu = [ciu,tgt];
    dff = sum(abs(bsxfun(@minus,ci,tgt))) == 0;
    mult = [mult sum(dff)];
    valu = [valu; max(val(dff))];
    indx{count} = c(dff);
    ci(:,dff) = [];
    val(dff) = [];
    c(dff) = [];
end
