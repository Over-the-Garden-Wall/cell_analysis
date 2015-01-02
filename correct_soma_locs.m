function correct_soma_locs(cell_nums)

    C = get_constants;
    close all
    
    for c = cell_nums
        
        soma_point = zeros(1,3);
        
        cell_dat = cell_data(c);
        current_mp = cell_dat.get_midpoint(true);
        
        fh = figure;
        plot_cells(c,2,.005,[0 0 0]); title([num2str(c) ': view 1']);
        hold on;
        scatter(current_mp(1), current_mp(3), '*r');
        
        [x1 y1] = getpts(fh);
        
        soma_point(3) = y1(end);
        
        close(fh);
        
        
        fh = figure;
        plot_cells(c,3,.005,[0 0 0]); title([num2str(c) ': view 2']);
        hold on;
        scatter(current_mp(1), current_mp(2), '*r');
        
        [x2 y2] = getpts(fh);
        
        soma_point(2) = y2(end);
        soma_point(1) = (x2(end)+x1(end))/2;
        
        
        close(fh);
        
        save([C.soma_dir 'cell_' num2str(c) '_soma.mat'], 'soma_point'); 
        
    end
    
    for c = cell_nums
        cell_dat = cell_data(c,true);
    end
        