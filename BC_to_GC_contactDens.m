function [bc_contact, bc_overlap, pr_pred] = BC_to_GC_contactDens(cell_nums, BC_nums, skip_factor, show_plot)

C = get_constants;


if ~exist('skip_factor', 'var') || isempty(skip_factor)
    skip_factor = 10;
end

if ~exist('show_plot', 'var') || isempty(show_plot)
    show_plot = false;
end

num_cells = length(cell_nums);


% BC_types = {'t1', 't2', 't3a', 't3b', 't4', 't5w', 't5l', 't5h', 't6', 't7'};
% type_layer = [1 1 1 1 1 2 2 2 2 2];

BCs = []; BC_typenum = []; 
for type_k = 1:length(BC_nums);
    BCs = [BCs BC_nums{type_k}];
    BC_typenum(end+1:length(BCs)) = type_k;
end
num_bc = length(BCs);

bc_contact = zeros(num_cells, num_bc);
bc_overlap = zeros(num_cells, num_bc);

for dk = 1:num_cells;
    d = cell_nums(dk);
        
    
    c_d = cell_data(d);
          
        d_contacts = double(c_d.contacts);
    
    p = c_d.get_surface;
    pdepth = C.f(p(:,1));
    
    p = p(pdepth < 100,2:3);
        
    p = p(ceil(1:skip_factor:end),:);
    
    tic
    for bk = 1:num_bc;
        b = BCs(bk);
        bc_data = cell_data(b);
        
        bc_hull = [];
        [bc_hull(:,1), bc_hull(:,2)] = poly2cw(bc_data.hull_2d(:,1), bc_data.hull_2d(:,2));
        
        
        bc_contact(dk, bk) = sum(d_contacts(2, d_contacts(1,:)==b));
        
        bc_overlap(dk, bk) = sum(inpolygon(p(:,1), p(:,2), ...
            bc_hull(:,1), bc_hull(:,2))) * skip_factor;
        
        
    end
    toc
    
end


pr_pred = peters_rule_prediction(cell_nums, BC_nums);

if show_plot

    densities = zeros(length(cell_nums), length(BC_nums));
    for n = 1:length(BC_nums)
        densities(:, n) = sum(bc_contact(:, BC_typenum==n),2) ./ sum(bc_overlap(:, BC_typenum==n),2);
    end
    
%     densities = bc_contact./bc_overlap;

    
    
    mean_dens = zeros(1,size(densities,2));
    std_dens = zeros(1,size(densities,2));

    for n = 1:size(densities,2)
        mean_dens(n) = mean(densities(~isnan(densities(:,n)),n));
        std_dens(n) = std(densities(~isnan(densities(:,n)),n));    
    end
    ste_dens = std_dens / sqrt(size(densities,1)-1);
    figure; 
    
    h = error_dot_plot(100*[mean_dens' pr_pred], 100*[ste_dens' zeros(length(pr_pred),1)]);

    legend(h, {'Observed Contact', 'Expected Contact'});
    prep_figure(gcf,gca)
    ylabel('Percent Coverage', 'fontSize', 12/.45);

end


end
