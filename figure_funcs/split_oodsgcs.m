C = get_constants;


oocns = C.type.oodsgc;
num_dsgc = length(oocns);
       
for dn = 1:num_dsgc
    
    d = oocns(dn);
    
    c_d = cell_data(d);
    
    
    d_p = c_d.get_surface;
    
    p_depth = C.f(d_p(:,1));
    
    p{1} = d_p(p_depth > 10 & p_depth < 50, 1:3);
    p{2} = d_p(p_depth > 50 & p_depth < 80, 1:3);
    
%     layer_mid = zeros(2,2);
    
    all_conts = double(c_d.contacts);
    
    
    hulls = cell(2, 1);

    
    for l = 1:2
        
% 
%         hull_inds = convhull(p{l}(:,1),p{l}(:,2));
%         hulls{l} = [p{l}(hull_inds,1) p{l}(hull_inds,2)];
        
        
        soma_point = mean(p{l});
        surface_points = p{l};
        
        save([C.point_dir 'cell_' num2str(d + l*1000) '_surface.mat'], 'surface_points'); 
        save([C.soma_dir 'cell_' num2str(d + l*1000) '_soma.mat'], 'soma_point');
        
        cell_data(d + l*1000, true);
    end
    
end