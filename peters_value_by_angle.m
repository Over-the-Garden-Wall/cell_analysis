function grid_bin_peters_script

    C = get_constants;
%     cell_types = {'t1', 't2', 't3a', 't3b', 't4'};
    cell_types = {'off_sac'};
    thresh = 10^7;
%     targ_nums = C.type.off_sac;
    targ_nums = 10010; %C.type.j;
    num_targs = length(targ_nums);
    
    cell_dats = cell(num_targs,1);    
    for n = 1:num_targs
        cell_dats{n} = cell_data(targ_nums(n));
    end
    
%     figure;
%     all_ax = gca;
%     hold all
%     
%     sub_f = figure;
%     sub_h = zeros(length(cell_types),1);
%     for t = 1:length(cell_types)
%         sub_h(t) = subplot(2,3,t);
%     end
%     
    all_ax = [];
    cols = colormap('Lines');
    h = zeros(length(cell_types),1);
    
    for t = 1:length(cell_types);
        
        cns = C.type.(cell_types{t});
        num_cells = length(cns);
        
        overlap = zeros(num_cells,num_targs);
        contact = zeros(num_cells,num_targs);
        angles = zeros(num_cells,num_targs);
        
        off_dats = cell(num_cells,1);
        for cn = 1:num_cells
            off_dats{cn} = cell_data(cns(cn));
        end
        
        for n = 1:num_targs
            targ_loc = cell_dats{n}.get_midpoint;
            for cn = 1:num_cells
                off_loc = off_dats{cn}.get_midpoint;
                overlap(cn,n) = compute_cell_overlap([targ_nums(n) cns(cn)]);

                if cell_dats{n}.contact_map.isKey(cns(cn));
                    cont_id = cell_dats{n}.contact_map(cns(cn));

                    contact(cn,n) = cell_dats{n}.contact_area(cont_id);
                end
                
                angles(cn,n) = atan2(off_loc(3)-targ_loc(3),off_loc(2)-targ_loc(2));


            end
        end
        
        contact = contact(overlap(:)>=thresh);
%         label_cns = cns(overlap(:)~=0);
        angles = angles(overlap(:)>=thresh);
        overlap = overlap(overlap(:)>=thresh);
        
        peters_ind = contact./overlap;
        
        
%          [corr_p] = corr(contact, overlap);
            
%             figure;
        if isempty(all_ax)
            figure;
            h(t) = polar_scatter(angles(:), peters_ind, '*', 'markerFaceColor', cols(t,:));
            all_ax = gca; hold all;
        else
            h(t) = polar_scatter(all_ax, angles(:), peters_ind, '*', 'markerFaceColor', cols(t,:));
        end
        
      
%         fit_slope = mean(contact./overlap);
        
        
%         linearCoef = polyfit(overlap,contact,1);
%         linearFit = polyval(linearCoef,overlap);
%         plot(overlap,linearFit,'-', 'Color', cols(t,:));
%         plot(all_ax, overlap,fit_slope*overlap,':', 'Color', cols(t,:));
        
        
%         set(sub_f, 'CurrentAxes', sub_h(t));
%         scatter(sub_h(t), overlap,contact,'*', 'markerFaceColor', cols(t,:));
%         hold(sub_h(t), 'on');
%         for cn = 1:length(overlap)
% %                 
%             text(overlap(cn), contact(cn) + max(contact)/30, num2str(label_cns(cn)-70000));
% % %             SA = get_size_stats(cell_nums(n));
% % %             all_cont(n) = all_cont(n)/SA;
%         end
        
        
%         title(sub_h(t), [cell_types{t} ', p=' num2str(corr_p)]);
    end
    legend(h, cell_types);
end