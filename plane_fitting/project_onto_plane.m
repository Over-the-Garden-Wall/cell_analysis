function [new_coords t] = project_onto_plane(coords, P)
    %coords is n x 3, P is 4x1

    r = sum(P(1:3).^2);
    t = -(coords*P(1:3)' + P(4))/r;
    
    %t is now 1 x n
    new_coords = coords + t*P(1:3);
end
    