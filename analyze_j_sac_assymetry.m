function analyze_j_sac_assymetry(j_num)

    C = get_constants;
%     load('/net/omicfs/home/matthew/stratification/j_synapses.mat');
    load(C.conn_loc)
    syn_conns = double([conns(2:6, conns(1,:)==j_num) conns([1 3:6], conns(2,:)==j_num)]);
    syn_conns = syn_conns(:,syn_conns(1,:)>=70000 & syn_conns(1,:)<80000);
    
    cell_dat = cell_data(j_num);
    
    p = cell_dat.get_surface;
    hull_p = convhull(p(:,2), p(:,3));
    
    hull = p(hull_p, 2:3);
    
    syn_hull_p = convhull(syn_conns(4,:)', syn_conns(5,:)');
    syn_hull = syn_conns(4:5,syn_hull_p)';
    
    
    
    figure;    
    plot_cells(j_num,1,.001,[0 0 0]); hold on
     plot(hull(:,1), hull(:,2), 'LineWidth', 2);
    hold all; plot(syn_hull(:,1), syn_hull(:,2), 'LineWidth', 2);
    
    
    num_syns = size(syn_conns,2);
    vec_begin = zeros(num_syns,2);
    vec_end = zeros(num_syns,2);
    syn_size = zeros(num_syns,1);
    syn_vec = zeros(num_syns,2);
    
    for s = 1:num_syns
        sac_dat = cell_data(syn_conns(1,s));
        sac_soma = sac_dat.get_midpoint(true);
        
        vec_end(s,:) = syn_conns(4:5,s)';
        
        [int_x int_y] = line_intersect(sac_soma(2:3), vec_end(s,:), hull(1:end-1,:), hull(2:end,:));
        if isempty(int_x)
            vec_begin(s,:) = sac_soma(2:3);
        else
            vec_begin(s,:) = [int_x int_y];
        end
        
        syn_vec(s,:) = (vec_end(s,:)-vec_begin(s,:));
        vec_r = sqrt(sum((syn_vec(s,:)).^2));
        syn_vec(s,:) = syn_vec(s,:)/vec_r;
        syn_size(s) = syn_conns(2,s)*vec_r;
        
        %       plot([hull_entry_loc(s,1); syn_conns(4,s)], [hull_entry_loc(s,2); syn_conns(5,s)], 'LineWidth', 2, 'Color', [0 0 0]);
%         scatter(sac_soma(2), sac_soma(3));
%         scatter(syn_conns(4,s), syn_conns(5,s));
%         
    end
    scatter(vec_end(:,1),vec_end(:,2), '*');
    
    j_mid = cell_dat.get_midpoint(true);
    mean_vec = sum([syn_vec(:,1).*syn_size, syn_vec(:,2).*syn_size])/10000;
    plot([j_mid(2); j_mid(2)-mean_vec(1)], [j_mid(3); j_mid(3)-mean_vec(2)], 'LineWidth', 2)
    
    
    
    rads = -pi:.01:pi;
    inhib = zeros(size(rads));
    angle_thresh = .6;
    for k = 1:length(rads);
        r = rads(k);
        r_vec = [sin(r), cos(r)];
        dprod = syn_vec(:,2)*r_vec(1) + syn_vec(:,1)*r_vec(2);
        dprod(dprod < angle_thresh) = 0;
        inhib(k) = sum(dprod.*syn_size);
    end
    figure;
    polar(rads, inhib);
    
end