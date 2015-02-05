C = get_constants;
oocns = C.type.oodsgc;
colmap = make_colormap(13);

for n = 1:length(oocns)
    cn = oocns(n);
    
    c_d = cell_data(cn);
    soma_loc = c_d.get_midpoint(true);
    figure; 
    set(gcf, 'Position', [0 0 600 1200]);
    
    subplot(1,2,1); hold on
    plot_cells(cn, 1, .01, colmap(n,:));
    scatter(soma_loc(2), soma_loc(3), '*', 'markerEdgeColor', [0 0 0]);
    
    subplot(1,2,2); hold on
    plot_cells(cn, 3, .01, colmap(n,:));
    
    set(gcf, 
    
    
    