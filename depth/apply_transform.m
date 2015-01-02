function p = apply_transform(T,p)

    p = p*T.global_Q';
    
    dist = sqrt((p(2)-T.nodes(:,1)).^2 + ...
        (p(3)-T.nodes(:,2)).^2);

    valid_inds = dist<= T.reach;
%     valid_inds = (dist==min(dist));
    
    if ~any(valid_inds)
        [~, d_ind] = min(dist);
        q = T.quaternion(d_ind,:);
        translatn = T.translation(d_ind,:);
    else
        weight = T.reach - dist(valid_inds);
        weight = weight/sum(weight);                

        q = sum(T.quaternion(valid_inds,:).*(weight*ones(1,4)),1);
        translatn = sum(T.translation(valid_inds,:).*(weight*ones(1,3)), 1);
    end

%     Q = angles2Q(theta, phi, psi);
    Q = quaternion2Q(q);

    p = p*Q' + translatn;
end