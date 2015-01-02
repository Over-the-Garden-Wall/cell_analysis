% types = {'t1', 't2', 't3a', 't3b', 't4'};
% 
% x_f = @(x,y)(x(ceil(size(x,1)*.75),1));
% y_f = @(x,y)(x(ceil(size(x,1)*.25),1));

% 
% 
% types = {'t1', 't2'};
% 
% x_f = @(x,y)(x(ceil(size(x,1)*.75),1));
% y_f = @(x,y)(x(ceil(size(x,1)*.25),1));



types = {'t3a', 't3b', 't4'};

x_f = @(x,y)(x(ceil(size(x,1)*.1),1));
% x_f = @(x,y)(y.hull_area);
y_f = @(x,y)(y.V);



num_types = length(types);

C = get_constants;

x_data = cell(num_types,1);
y_data = cell(num_types,1);

figure; hold on

for k = 1:length(types)
    num_cells = length(C.type.(types{k}));
    x_data{k} = zeros(num_cells,1);
    y_data{k} = zeros(num_cells,1);
    
    
    for ck = 1:num_cells
        c = C.type.(types{k})(ck);
        cell_dat = cell_data(c);
        p = cell_dat.get_surface;
        d = sort(C.f(p(:,1)));
        d(d<0) = [];

        x_data{k}(ck) = x_f(d,cell_dat);
        y_data{k}(ck) = y_f(d,cell_dat);
    end

    scatter(x_data{k}, y_data{k}, '*', 'MarkerEdgeColor', C.colormap(k,:));
    
    
end

prep_figure(gcf,gca)