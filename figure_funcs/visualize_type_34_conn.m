function visualize_type_34_conn()

    C = get_constants;

    cell_nums = [C.type.t3a, C.type.t3b, C.type.t4];
    num_cells = length(cell_nums);   
    
    contact = zeros(num_cells,1);
    perc10 = zeros(num_cells,1);
    vol = zeros(num_cells,1);
    
    
    for n = 1:num_cells;
        cell_dat = cell_data(cell_nums(n));
        
        vol(n) = cell_dat.V;
        
        p = cell_dat.get_surface;
        d = C.f(p(:,1));
        d = sort(d);
        d(d<0) = [];
        perc10(n) = d(ceil(.1*length(d)));
        
        cont = double(cell_dat.contacts);
        cont = cont(:,cont(1,:)>=70000 & cont(1,:) < 80000);
        num_sacs = length(unique(cont(1,:)));
        
        contact(n) = sum(cont(2,:))/num_sacs;
    end
    
%     log_cont = log(1+contact);
    figure; hold on;
    for n = 1:num_cells
        scatter(perc10(n), vol(n), '*', 'MarkerEdgeColor', .8*([1 1 1] - contact(n)/max(contact)));
    end
end
        