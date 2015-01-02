function [Q phi psi] = find_planar_rotation(P)
    
    if P(3) == 0
        phi = 0;
        Q = eye(3);
    else
        phi = -atan2(P(3),P(1));
        Q = [cos(phi) 0 -sin(phi); ...
            0 1 0; ... 
            sin(phi) 0 cos(phi)];
    end
    
    P = Q*P(1:3)';
    
    if P(2) ~= 0
        psi = -atan2(P(2),P(1));
        Q = [cos(psi) -sin(psi) 0; ...
            sin(psi) cos(psi) 0; ...
            0 0 1]*Q;
    else
        psi = 0;
    end
    
    
end

        
    