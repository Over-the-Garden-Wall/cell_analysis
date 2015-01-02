function [t f s] = Q2angles(Q)
     
%     Q = [ct*cs, sf*st*cs-cf*ss, sf*ss+cf*st*cs; ...
%         ct*ss, cf*cs + sf*st*ss, cf*st*ss-sf*cs; ...
%         -st, sf*ct, cf*ct];
    acceptable_error = .00001;


    t = asin(-Q(3,1));
    f = asin(Q(3,2)/cos(t));
    s = acos(Q(1,1)/cos(t));
    
    %resolve sign differences
    
    %assume t is correct
    if abs(Q(3,3) - cos(t)*cos(f)) > acceptable_error
        f = pi - f;
    end
    if abs(Q(2,1) - cos(t)*sin(s)) > acceptable_error
        s = -s;
    end    
    %will give correct Q... unless t is wrong
    
    if abs(Q(2,2) - cos(f)*cos(s) - sin(f)*sin(t)*sin(s)) > acceptable_error
        t = pi-t;
        f = asin(Q(3,2)/cos(t));
        s = acos(Q(1,1)/cos(t));
        if abs(Q(3,3) - cos(t)*cos(f)) > acceptable_error
            f = pi - f;
        end
        if abs(Q(2,1) - cos(t)*sin(s)) > acceptable_error
            s = -s;
        end    
    end
        
    
    
end