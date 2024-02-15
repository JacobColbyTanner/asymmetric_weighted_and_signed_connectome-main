function [cinew,idx] = fcn_sort_communities(ci)
%FCN_SORT_COMMUNITIES   reorder community labels by size
%
%   [cinew,idx] = fcn_sort_communities(ci) takes a partition, ci, and
%   relabels its communities so that the largest community is labeled 1,
%   next largest labeled 2, next largest 3, and so on.
%
%   Inputs:
%               ci,     [node x 1] partition
%
%   Outputs:
%
%            cinew,     [node x 1] relabeled partition
%              idx,     communities in the original partition ordered by
%                       size (descending order)
%
%   Richard Betzel, Indiana University, 2020
%   
h = hist(ci,1:max(ci));
[~,idx] = sort(h,'descend');
cinew = zeros(size(ci));
for j = 1:max(ci)
    jdx = ci == idx(j);
    cinew(jdx) = j;
end