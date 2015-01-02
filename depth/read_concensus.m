function segs = read_concensus(fn)

    fid = fopen(fn);

    M = textscan(fid, '%u32 %*[^\n]');
    
    fclose(fid);
    
    segs = M{1};

end