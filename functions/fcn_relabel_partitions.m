function cinew = fcn_relabel_partitions(ci)
%FCN_RELABEL_PARTITIONS     relabels partitions
%
%   cinew = fcn_relabel_partitions(ci) takes as input an ensemble of
%   partitions, ci, of dimensions [nodes x partitions] and relabels
%   communities in the order in which each community is encountered. Note
%   that the relabeling is arbitrary, but it facilitates
%   identifying repeats of the same partition and can sometimes assist with
%   visualization.
%
%   Inputs:
%           ci,     [node x partition] ensemble of partitions
%
%   Outputs:
%        cinew,     [node x partition] ensemble or relabeled partitions.
%
%   Richard Betzel, Indiana University, 2012

[n,m] = size(ci);

cinew = zeros(n,m);
for i = 1:m
    c = ci(:,i);
    d = zeros(size(c));
    count = 0;
    while sum(d ~= 0) < n
        count = count + 1;
        ind = find(c,1,'first');
        tgt = c(ind);
        rep = c == tgt;
        d(rep) = count;
        c(rep) = 0;
    end
    cinew(:,i) = d;
end