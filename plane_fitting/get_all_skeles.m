function data = get_all_skeles()

    on_cells = [4,5,6,10,11,12,15,16,17,27,33];

    dirs = dir('./');
    
    k = 0;
    for n = 1:length(dirs)
        if length(dirs(n).name)>2 && strcmp(dirs(n).name(1:3),'sac')
            k=k+1;
        end
    end
    
    data = cell(k,1);
    k = 0;
    for n = 1:length(dirs)
        if length(dirs(n).name)>2 && strcmp(dirs(n).name(1:3),'sac')
            k=k+1;
            data{k} = analyze_skeleton(dirs(n).name);
            data{k}.sac_num = str2double(dirs(n).name(4:end));
            if any(on_cells==data{k}.sac_num)
                data{k}.on_cell = true;
            else
                data{k}.on_cell = false;
            end
        end
    end
end