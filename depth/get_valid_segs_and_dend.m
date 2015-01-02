function uni_sprvox = get_valid_segs_and_dend(cell_id, out_dir)

    if ~exist(['./' out_dir],'dir')
        mkdir(['./' out_dir]);
    end

    disp(['working on ' cell_id]);
    tic
    
    [seg_list_fns, omni_file_fns] = get_filenames_for_cell(cell_id);
    
    num_files = length(seg_list_fns);
    
    dends = cell(num_files,1);
    dendVals = cell(num_files,1);
    valid_segs = cell(num_files,1);
    uni_sprvox = cell(num_files,1);
    point_count = cell(num_files,1);
    
    for n = 1:num_files
        valid_segs{n} = read_concensus(seg_list_fns{n});
        [dends{n} dendVals{n}] = read_omni_dendrogram(omni_file_fns{n});
        uni_sprvox{n} = unique([dends{n}(:); valid_segs{n}(:)]);
        
        point_count{n} = count_points_omni(omni_file_fns{n}, uni_sprvox{n});
        
    end
    
    
    save(['./' out_dir '/' cell_id '_info.mat'],'valid_segs', 'dendVals', 'dends', 'uni_sprvox', 'point_count', 'omni_file_fns');
    
    toc
end
    