function [start_coords end_coords] = omni_fn2coords(fn)

    %example: e2198_z_s16_131_91_e20_146_99.omni.files
    
    s_loc = strfind(fn,'s');
    fn = fn(s_loc(1)+1:end);
    per_loc = strfind(fn,'.');
    fn = fn(1:per_loc(1)-1);
    
    e_loc = strfind(fn,'e');
    
    start_substring = fn(1:e_loc(1)-2);
    us_loc = strfind(start_substring, '_');
    
    start_sectors = [str2double(start_substring(1:us_loc(1)-1)), ...
        str2double(start_substring(us_loc(1)+1:us_loc(2)-1)), ...
        str2double(start_substring(us_loc(2)+1:end))];
    
    end_substring = fn(e_loc+1:end);
    us_loc = strfind(end_substring, '_');
    
    end_sectors = [str2double(end_substring(1:us_loc(1)-1)), ...
        str2double(end_substring(us_loc(1)+1:us_loc(2)-1)), ...
        str2double(end_substring(us_loc(2)+1:end))];
    
    
    %128*(X-1)+18    128*Y+18,    128*Z+18
    
    start_coords = start_sectors*128 + 18 - [128 0 0];
    end_coords = end_sectors*128 + 17 + [0 128 128];
    
    
    
end
    
    
    
        