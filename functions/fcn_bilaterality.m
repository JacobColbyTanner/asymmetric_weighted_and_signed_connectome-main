function lat = fcn_bilaterality(c,hemi,nrand)
%FCN_BILATERALITY    measures laterality of partition
%
%   lat = fcn_bilaterality(c,hemi,nrand) takes a [node x 1] partition, a
%   [node x 1] vector of hemisphere labels (1 = RH, 2 = LH), and a number
%   of random permutations. It returns the laterality index from Lohse et
%   al (2013, PLOS Comp Biol). Values of 0 correspond to balance --
%   communities are composed of equal numbers of nodes from both 
%   hemispheres. Values of 1 correspond to imbalance -- communities tend to
%   contain mostly nodes from one or the other hemisphere.
%
%   Inputs:
%             c,     [node x 1] vector of community labels
%          hemi,     [node x 1] vector of hemisphere labels
%         nrand,     scalar value for number of random permutations
%
%   Outputs:
%     
%           lat,     laterality index
%           
%   Example:
%       nrand = 1000;                           % number of randomizations
%       hemi = ones(size(sc,1),1);              % hemisphere labels
%       hemi(201:400) = 2;
%       gamma = 1;                              % resolution parameter
%       [c,q] = community_louvain(sc);          % estimate communities
%       lat = fcn_bilaterality(c,hemi,nrand);   % calculate laterality
%
%
%   Reference:
%       Lohse, C., Bassett, D. S., Lim, K. O., & Carlson, J. M. (2014). 
%       Resolving anatomical and functional structure in human brain 
%       organization: identifying mesoscale organization in weighted 
%       network representations. PLoS computational biology, 10(10), 
%       e1003712.
%
%   Jacob Tanner, Richard Betzel, Indiana University, 2024

%% calculate laterality

% create dummy variables for hemisphere ...
H = dummyvar(hemi);

% ... and community vectors
C = dummyvar(c);

% calculate contingency table
P = H'*C;

% size of each community
sz = sum(C);

% calculate uncorrected laterality as |N_LH - N_RH|/|N_C| for each
% community
L = abs(diff(P)./sz);

% calculate expected laterality -- note that the random permutation could
% be replaced with geometry preserving models, e.g. spin tests.
n = length(c);
Lr = zeros(max(c),nrand);
for irand = 1:nrand
    r = randperm(n);
    Hr = H(r,:);
    Pr = Hr'*C;
    Lr(:,irand) = abs(diff(Pr)./sz);
end

% scale community laterality by size
lat = L.*sz;

% do the same for null model
latnull = mean(bsxfun(@times,Lr,sz'),2);

% calculate difference and divide by number of nodes
lat = (sum(lat) - sum(latnull))/n;