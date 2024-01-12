function b = get_geometric_null(sc)

    
    mtrx = sc;
    mask = mtrx ~= 0;
    g = nanmean(mtrx(mask));
    b = mtrx - (mask*g);