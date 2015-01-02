function im = thin2full(cell_no, dsmp_fact)
   
    root_dir = '~/stratification/thinned_data/thin_'; 

    fn = [root_dir num2str(cell_no) '.mat'];
    
    
    load(fn);
    
    max_size = max(p);
%     dsmp_fact = uint32(dsmp_fact);
    
    max_size = floor(max_size./dsmp_fact+1);
    
    im = false(max_size);
    
    for k = 1:3;
        p(:,k) = floor(p(:,k)/dsmp_fact(k)+1);
    end
    
    
        
        im(sub2ind(max_size,p(:,1),p(:,2),p(:,3))) = true;
    
    
end
    
    

    