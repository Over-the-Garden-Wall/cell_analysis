function run_length_approx1(cell_type_char, cell_type_num)

    load(['./forSrini/' cell_type_char '_info.mat'])

    total_lengths = get_lengths_for_cell_blocks(cell_type_num, cell_type_char);
    
    figure; 
    
    c = colormap('Lines');
    
    for n = 1:length(dends)
        


        [total_pieces change_type] = analyze_dends(dends{n}, dendVals{n}, valid_segs{n}, point_count{n}, uni_sprvox{n});
            
        
        
        false_positives = cumsum(change_type(:,2));
        run_length = cumsum(change_type(:,1))/total_pieces*total_lengths;
            
        
        plot(false_positives, run_length, 'Color', c(n,:));
            hold on;
        
    end
    
end
        