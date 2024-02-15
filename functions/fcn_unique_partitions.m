function [ciu,valu,mult,indx] = fcn_unique_partitions(ci,val)
%FCN_UNIQUE_PARTITIONS      returns unique partitions
%
%   [ciu,valu,mult,indx] = fcn_unique_partitions(ci,val) takes a set of
%   partitions ci and associated scalar values val and identifies only the
%   unique partitions, ciu, the max val for each ciu, and the multiplicity
%   mult of that partition.
%
%   Inputs:     ci,     partitions
%               val,    scalar values
%
%   Outputs:    ciu,    unique partitions
%               valu,   associated values
%               mult,   multiplicity
%               indx,   indices
%
%   Richard Betzel, Indiana University, 2012
%

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