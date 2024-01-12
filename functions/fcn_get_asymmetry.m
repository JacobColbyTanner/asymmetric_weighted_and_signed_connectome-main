function [asymmetry, abs_asymmetry] = fcn_get_asymmetry(B)


    asymmetry = B-B';
    abs_asymmetry = abs(B-B');

