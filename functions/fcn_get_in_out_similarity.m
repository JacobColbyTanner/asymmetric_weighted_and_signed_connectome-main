function [inout_sim,p] = fcn_get_in_out_similarity(B)

    id = B == 0;
    B(id) = nan;
    
    for i = 1:size(B,1)
        
        out = B(i,:)';
        inn = B(:,i);

        [r(i),p(i)] = corr(out,inn,"rows","complete");

    end
    

    inout_sim = r;


