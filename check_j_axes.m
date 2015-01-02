function check_j_axes
    close all
    C = get_constants;
    k = 1;
    for n = C.type.j
        cell_dat = cell_data(n);
        
        subplot(1,2,k);
        plot_cells(n,1,.001, [0 0 0]);
        hold all
        
        soma = cell_dat.get_midpoint(true);
        
        dist_point = soma(2:3) + 100000*cell_dat.dist_axis;
        
        plot([soma(2); dist_point(1)], [soma(3); dist_point(2)],'Color', [1 0 0],'LineWidth', 3);
        
        soma = round(soma); dist_point = round(dist_point);
        disp([num2str(n) '    ' num2str(soma(1)) '    ' num2str(soma(2)) '    ' num2str(soma(3)) '    ' num2str(dist_point(1)) '    ' num2str(dist_point(2))])
        
        title(num2str(n));
        k = k+1;
    end
end