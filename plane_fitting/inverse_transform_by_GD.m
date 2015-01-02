function p = inverse_transform_by_GD(desired_p, T)
    C = get_constants;
    eta = .00001;
    delta = .0001;
    p = [500 500 500];
    
    err_func = @(x,y) sum(abs(apply_transform(T,x)-y).^2);
    
    E = err_func(p, desired_p);
%     new_E = zeros(3,1);
    
    dEdp = zeros(1,3);
    
    while E>2
    
        for d = 1:3
            p(d) = p(d)+delta;
            new_E = err_func(p, desired_p);
%           new_E = sum(abs(apply_transform(T,p)-desired_p).^1.1);\
            p(d) = p(d)-delta;
            dEdp(d) = (new_E-E)/delta;
            
        end
            
        
        p = p-eta*dEdp;
        
        
        E = err_func(p, desired_p);
%         disp(E)
    end
end
        
        