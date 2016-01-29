C = get_constants;

sacs = 'on_sac';
hull_cells = {'BC7', 'BC5t'};

hulls = cell(length(hull_cells),2);

for n = 1:length(hull_cells)
    for c = C.type.(hull_cells{n})
        c_d = cell_data(c);
        h = [];
        [h(:,1), h(:,2)] = poly2cw(c_d.hull_2d(:,1), c_d.hull_2d(:,2));
        [hulls{n,:}] = polybool('union', hulls{n,1}, hulls{n,2}, h(:,1), h(:,2));
    end
end

for n = length(hull_cells):-1:2
    [hulls{1,:}] = polybool('intersection', hulls{n,1}, hulls{n,2}, hulls{1,1}, hulls{1,2});
end

cns = [];
for c = C.type.(sacs)
    new_cell_nums = subarbor_cell_creation(c, [hulls{1,1} hulls{1,2}], 1000);
    cns = [cns new_cell_nums];
end

fprintf('\n');
fprintf('%d ', cns);
fprintf('\n');