function [edge_usage,percent_usage] = fcn_get_edge_usage(sc)

has_neg = sum(sc<0,"all") > 0;
directed = sum(sc-sc'~=0,"all") > 0;

if has_neg == 0
    L = sparse(1./sc);
elseif has_neg == 1
    sf = 0.0001;
    epsilon = range(nonzeros(sc))*sf;
    L = sc;
    L(sc ~= 0) = L(sc ~= 0) - min(sc(:));% + epsilon;
    L(sc ~= 0) = L(sc ~= 0) + epsilon;
    L = 1./L;
end


n = length(L);
if directed == 1
    G = digraph(L);
elseif directed == 0
    G = graph(L);
end

edge_usage = zeros(n);
for source = 1:n
    [pth_all] = shortestpathtree(G,source,'OutputForm','cell');
    

    for target = 1:n
        
        if source ~= target
            pth = pth_all{target};
            for j = 1:length(pth) - 1
                edge_usage(pth(j),pth(j + 1)) = edge_usage(pth(j),pth(j + 1)) + 1;
            end

        end
    end
end

total_edges = sum(sc~=0,"all");
used = sum(edge_usage > 0,"all");

percent_usage = used/total_edges;


