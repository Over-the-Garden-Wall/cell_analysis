function Q = get_best_rigid_transform()
    data = get_all_skeles;
    
    
    Ps = zeros(length(data),4);
    weight = zeros(length(data),1);
    
    for n = 1:length(data);
        Ps(n,:) = data{n}.P(1:4);
        weight(n) = size(data{n}.coords,1);        
    end
    
    P = sum(Ps .* (weight*ones(1,4)))/sum(weight);
    P(1:3) = P(1:3)/sqrt(sum(P(1:3).^2));
    
    Q = find_planar_rotation(P);
    
end
    
    