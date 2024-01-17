function [asymmetry, abs_asymmetry] = fcn_get_asymmetry(B)

%Function to measure the laterality of a modular partition. Laterality measures how much a module favors one hemisphere over the other.
%
%   
%
%   Inputs:
%       B,
%           modular partition of each node into different communities represented as a list of community labels
%       
%              
%
%   Outputs:
%     
%       asymmetry,
%           asymmetry matrix where each value gives the difference between edge(i,j) and edge(j,i)
%       abs_asymmetry,
%           absolute asymmetry matrix where each value gives the absolute difference between edge(i,j) and edge(j,i)
%           
%
%   Example:
%      
%        [asymmetry, abs_asymmetry] = fcn_get_asymmetry(B)


    asymmetry = B-B';
    abs_asymmetry = abs(B-B');

