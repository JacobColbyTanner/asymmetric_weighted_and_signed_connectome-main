function lat = fcn_bilaterality(c,hemi,nrand)
% 
% clear all
% close all
% clc
% 
% load ../mat/bmat.lag=1.samplegamma.N_92.mat
% c = cisample(:,1);
% hemi = [ones(46,1); ones(46,1)*2];
% nrand = 1000;


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