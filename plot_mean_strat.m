function strats = plot_mean_strat(cell_nos, relevant_portion, invertaxis, colmap)



    if ~exist('invertaxis','var') || isempty(invertaxis);
        invertaxis = 'none';
    end
    if ~exist('colmap', 'var') || isempty(invertaxis)
        colmap = colormap('Lines');
    end
    
    C = get_constants;
%     colmap = C.colormap;
%     colmap = make_colormap(length(cell_nos), 5);
    
    
    if ~exist('relevant_portion','var');
        relevant_portion = [C.strat_x(1) C.strat_x(end)];
    end

    num_types = length(cell_nos);
    
    strats = zeros(length(C.strat_x),num_types);
    
    for k = 1:num_types
    
        for cell_no = cell_nos{k}

            c_d = cell_data(cell_no);

            C = get_constants;

            s = c_d.stratification;
            if length(s) > size(strats,1)
                s = s(1:size(strats,1));
            end
%             disp([length(s) size(strats,1)]);
            strats(1:length(s),k) = strats(1:length(s),k) + s/length(cell_nos{k});
            
        end
        
    end
    
    strats = strats((relevant_portion(1) - C.strat_x(1) + 1):(relevant_portion(2) - C.strat_x(1) + 1),:);
    for k = 1:num_types
        strats(:,k) = strats(:,k)/sum(strats(:,k));
    end
    
    figure; hold all;
     if strcmp(invertaxis,'none')
         for k = 1:num_types   
         plot(relevant_portion(1):relevant_portion(end), strats(:,k), 'lineWidth', 2, 'Color', colmap(k,:));
         end
         set(gca, 'YTick', []);
     elseif strcmp(invertaxis,'xy')
            plot(strats(:,k), relevant_portion(1):relevant_portion(end), 'lineWidth', 2);
            set(gca, 'Ydir', 'reverse');
            
            figure;
            plot(cumsum(strats(:,k))*100, relevant_portion(1):relevant_portion(end), 'lineWidth', 2);
            set(gca, 'Ydir', 'reverse');
            set(gca, 'XTick', []);
     elseif strcmp(invertaxis,'x')
         for k = 1:num_types   
         plot(relevant_portion(1):relevant_portion(end), strats(:,k), 'lineWidth', 2);
         end
         set(gca, 'Xdir', 'reverse');
         set(gca, 'YTick', []);
        end
   
end