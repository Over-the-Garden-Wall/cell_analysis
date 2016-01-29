function h = ellipse_pca(p)

    mid_p = mean(p);
    p = p - ones(size(p,1),1)*mid_p;

    C = p'*p;
    
    [vecs, vals] = svd(C);
    
    angle = atan2(vecs(1,2), vecs(1,1));
    
    h = ellipse(sqrt(vals(1,1)/size(p,1))*2, sqrt(vals(2,2)/size(p,1))*2, angle, mid_p(1), mid_p(2));
end