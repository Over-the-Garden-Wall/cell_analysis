C = get_constants;


figure; ax_h = gca; plot_cells(10010,1,.001, .7*ones(1,3)); hold on

types = {'t1', 't2', 't3a', 't3b', 't4'};
% types = {'t1', 't4'};
num_types = length(types);


pos = cell(num_types,1);

c = colormap('Lines');

bins = 5:10:95;
j_dat = cell_data(10010);
soma = j_dat.get_midpoint(true);
soma = soma(2:3);
ax = j_dat.dist_axis;
perp_ax = [-ax(1), ax(2)];

cell_length = 117850;

for n = 1:length(bins)
    line_point = soma + ax*cell_length/100*bins(n);
    line_point1 = line_point([2 1]) - perp_ax*50000;
    line_point2 = line_point([2 1]) + perp_ax*50000;
    plot([line_point1(2); line_point2(2)], [line_point1(1); line_point2(1)], '--k');
end


for n = 1:num_types;
    num_cells = length(C.type.(types{n}));
    pos{n} = zeros(num_cells,3);
    name_mat = cell(num_cells,1);
    for k = 1:num_cells
        cell_dat = cell_data(C.type.(types{n})(k));
        pos{n}(k,:) = cell_dat.get_midpoint(false);
        
        
        
        if j_dat.contact_map.isKey(C.type.(types{n})(k))
            cont_area = j_dat.contact_area(j_dat.contact_map(C.type.(types{n})(k)));
        else
            cont_area = 0;        
        end
        name_mat{k} = [num2str(C.type.(types{n})(k)) ' ' num2str(cont_area)];
    end
    scat_hand(n) = scatter(pos{n}(:,2), pos{n}(:,3), '*', 'MarkerFaceColor', c(n,:));
%     gname(name_mat{k});
end
%     gname(name_mat);
        

% xlim = get(ax_h,'XLim');
% ylim = get(ax_h,'YLim');

    

legend(scat_hand, types);