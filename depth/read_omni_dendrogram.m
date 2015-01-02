function [dend, dendVals] = read_omni_dendrogram(omnifn)
    mst_file = [omnifn '/segmentations/segmentation1/segments/mst.data'];
    
    if ~exist(mst_file,'file')
        mst_file = [omnifn '/users/_default/segmentations/segmentation1/segments/mst.data'];
    end
    
    fid = fopen(mst_file);
    int_vals = fread(fid,'uint32');
    fclose(fid);
    
    fid = fopen(mst_file);
    dbl_vals = fread(fid,'double');
    fclose(fid);
    
    
    num_entries = length(dbl_vals)/4;
    
    dend = zeros(num_entries,2,'uint32');
    dendVals = zeros(num_entries,1);
%     
%     disp(num_entries)
%     disp(length(dbl_vals))
%     disp(length(int_vals))
    
    for n = 1:num_entries; 
        dend(n,:) = int_vals((n-1)*8+(2:3)); 
    end
    
    for n = 1:num_entries; 
        dendVals(n) = dbl_vals(4*(n-1)+3); 
    end
    
end