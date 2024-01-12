function [sign_asym] = fcn_get_sign_asymmetry(B)


    pos = B > 0;
    neg = B < 0;
    %get elements for i,j is pos/neg and j,i is neg,pos
    sign_asym = pos+neg' == 2;
