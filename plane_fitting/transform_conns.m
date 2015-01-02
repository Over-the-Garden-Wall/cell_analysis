function conns = transform_conns(conns,T)

    if size(conns,1) > 6
        conns(3,:) = sum(conns(3:5,:));
        conns(4:5,:) = [];
    end
    
    for n = 1:size(conns,2);
        conns(4:6,n) = int32(apply_transform(T,double(conns(4:6,n)'))');
    end
end