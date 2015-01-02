function P = fit_plane(coords)
    %c24*P = -c1

    coords(:,4) = 1;
    
    P = coords(:,2:4)\(-coords(:,1));
    
    P = [1 P'];    
end