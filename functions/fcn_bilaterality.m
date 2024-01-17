function lat = fcn_bilaterality(c,hemi,nrand)

   %Function to measure the laterality of a modular partition. Laterality measures how much a module favors one hemisphere over the other.
%
%   
%
%   Inputs:
%       c,
%           modular partition of each node into different communities represented as a list of community labels
%       hemi,
%           list of labels indicating which hemisphere each node belongs to
%       nrand,
%           number of randomly permuted null models to measure laterality against
%              
%
%   Outputs:
%     
%       lat,
%           laterality
%           
%
%   Example:
%       nrand = 1000;
%       hemi = ones(size(sc,1),1);
%       hemi(201:400) = 2: %label which nodes belong to each hemisphere
%       [c,Q] = community_louvain(sc,gamma,[],'neg_asym');
%       lat = fcn_bilaterality(c,hemi,nrand);
%       
%
%   


H = dummyvar(hemi);
C = dummyvar(c);
P = H'*C;
sz = sum(C);
L = abs(diff(P)./sz);


n = length(c);
Lr = zeros(max(c),nrand);
for irand = 1:nrand
    r = randperm(n);
    Hr = H(r,:);
    Pr = Hr'*C;
    Lr(:,irand) = abs(diff(Pr)./sz);
    
%     for i = 1:max(c)
%         h = hist(hemi(cr == i),1:2);
%         lambdacr(i,irand) = abs(h(2) - h(1))/sz(i);
%     end
end

lat = L.*sz;
latnull = mean(bsxfun(@times,Lr,sz'),2);

lat = (sum(lat) - sum(latnull))/n;

% lat = sum(lambdac'.*sz);
% latnull = mean(bsxfun(@times,lambdacr',sz),1);
% 
% lat = (sum(lat) - sum(latnull))/n;
