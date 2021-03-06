function grid_bin_peters_script

    C = get_constants;
    cell_types = {'t1', 't2', 't3a', 't3b', 't4'};
%     cell_types = {'off_sac'};
    
    targ_nums = C.type.off_sac;
%     targ_nums = 10010; %C.type.j;
    num_targs = length(targ_nums);
    
    cell_dats = cell(num_targs,1);
    
    for n = 1:num_targs
        cell_dats{n} = cell_data(targ_nums(n));
    end
    
    figure;
    all_ax = gca;
    hold all
    
    sub_f = figure;
    sub_h = zeros(length(cell_types),1);
    for t = 1:length(cell_types)
        sub_h(t) = subplot(2,3,t);
    end
    
    
    cols = colormap('Lines');
    h = zeros(length(cell_types),1);
    
    for t = 1:length(cell_types);
        
        cns = C.type.(cell_types{t});
        num_cells = length(cns);
        
        overlap = zeros(num_cells,num_targs);
        contact = zeros(num_cells,num_targs);
        
        name_mat = [];
        
        for n = 1:num_targs
            for cn = 1:length(cns)

                overlap(cn,n) = compute_cell_overlap([targ_nums(n) cns(cn)]);

                if cell_dats{n}.contact_map.isKey(cns(cn));
                    cont_id = cell_dats{n}.contact_map(cns(cn));

                    contact(cn,n) = cell_dats{n}.contact_area(cont_id);
                end

                name_mat{end+1} = num2str(cns(cn));
            end
        end
        
        contact = contact(overlap(:)~=0);
%         label_cns = cns(overlap(:)~=0);
        overlap = overlap(overlap(:)~=0);
        
        
         [corr_p] = corr(contact, overlap);
            
%             figure;
      h(t) = scatter(all_ax, overlap,contact,'*', 'markerFaceColor', cols(t,:));
        
      
        fit_slope = mean(contact./overlap);
        
        
%         linearCoef = polyfit(overlap,contact,1);
%         linearFit = polyval(linearCoef,overlap);
%         plot(overlap,linearFit,'-', 'Color', cols(t,:));
        plot(all_ax, overlap,fit_slope*overlap,':', 'Color', cols(t,:));
        
        
        set(sub_f, 'CurrentAxes', sub_h(t));
        scatter(sub_h(t), overlap,contact,'*', 'markerFaceColor', cols(t,:));
%         hold(sub_h(t), 'on');
%         for cn = 1:length(overlap)
% %                 
%             text(overlap(cn), contact(cn) + max(contact)/30, num2str(label_cns(cn)-70000));
% % %             SA = get_size_stats(cell_nums(n));
% % %             all_cont(n) = all_cont(n)/SA;
%         end
        
        
        title(sub_h(t), [cell_types{t} ', p=' num2str(corr_p) ', slope=' num2str(fit_slope)]);
    end
    legend(h, cell_types);
%     gname(name_mat);
end