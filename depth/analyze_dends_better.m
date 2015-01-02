function [total_pieces change_type change_val] = analyze_dends(dend_fn)

    load(dend_fn); %now have cell arrays: dendVals, dends, valid_segs
    
    num_vols = length(valid_segs);
    
    change_type = [];
    change_val = [];
    
    total_pieces = 0;
    for n = 1:num_vols
        total_pieces = total_pieces + length(valid_segs{n});
%         num_entries = size(dends,1);
        
%         introduces_error = false(num_entries,1);
        
        is_valid = sparse(double(valid_segs{n}),ones(length(valid_segs{n}),1),ones(length(valid_segs{n}),1), double(max(dends{n}(:))),1); 
        dend_is_in = full(is_valid(dends{n}));
        
        introduces_correct = dend_is_in(:,1) & dend_is_in(:,2);
        introduces_error = (dend_is_in(:,1) & ~dend_is_in(:,2)) | (dend_is_in(:,2) & ~dend_is_in(:,1));
        
        is_interesting = introduces_correct | introduces_error;
        change_val = [change_val; dendVals{n}(is_interesting)];
        change_type = [change_type; [introduces_correct(is_interesting) introduces_error(is_interesting)]];
    end
    
    [change_val, sort_inds] = sort(change_val, 'descend');
    change_type = change_type(sort_inds,:);
end
    
        