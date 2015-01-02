function spreadsheet_SAC_strats
    C = get_constants;
    bins = C.sac_bins;
    num_bins = length(bins);
    
    SACs = C.type.sure_off_sac;
    
    M = zeros(num_bins*length(SACs),102);
    for s = 1:length(SACs);
        cell_dat = cell_data(SACs(s));
        p = cell_dat.get_surface;
        d = C.f(p(:,1));
        sm = cell_dat.get_midpoint(true);
        
        r = sqrt((sm(2)-p(:,2)).^2+(sm(3)-p(:,3)).^2)/1.15;
        is_valid = d>0 & d<=99;
        r(~is_valid) = [];
        d(~is_valid) = [];
        
        bin_num = ceil(r/15000);
        
        for b = 1:num_bins
            k = (s-1)*num_bins+b;
            
            M(k,1) = SACs(s);
            M(k,2) = (b-.5)*15;
            
            M(k,3:102) = hist(d(bin_num==b),0:99);
        end
    end
    
    dlmwrite('./spreadsheets/SAC_strats_by_dist.txt',M);
    
end