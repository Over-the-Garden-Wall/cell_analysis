fns = dir('./point_data/');

cell_nums = zeros(length(fns),1);

for n = 1:length(fns)
    
    unders = strfind(fns(n).name,'_');
    if length(unders)<2
        cell_nums(n) = 0;
    else
        cell_nums(n) = str2double(fns(n).name(unders(1)+1:unders(2)-1));        
    end
    
end

stupid_M = cell(max(cell_nums),1);

for n = 1:length(fns)
    
    unders = strfind(fns(n).name,'_');
    if length(unders)<2
%         cell_nums(n) = 0;
    else
        cell_num = str2double(fns(n).name(unders(1)+1:unders(2)-1));  
        if cell_num > 0
        stupid_M{cell_num}{end+1} = fns(n).name;
        end
    end
    
end

for n = length(stupid_M):-1:1
    if isempty(stupid_M{n});
        stupid_M(n) = [];
    end
end

to_delete = {};

for n = 1:length(stupid_M)
    
    f_nums = zeros(length(stupid_M{n}),1);
    for k = 1:length(stupid_M{n});
        unders = strfind(stupid_M{n}{k},'_');
        f_nums(k) = str2double(stupid_M{n}{k}(unders(end)+1:end-4)); 
    end
    
    [f_nums sort_ind] = sort(f_nums);
    for k = 2:length(f_nums)
        if f_nums(k) ~= f_nums(k-1) + 1
            to_delete{end+1} = stupid_M{n}{sort_ind(k)};
        end
    end
    
end

for n = 1:length(to_delete)
    delete(['./point_data/' to_delete{n}]);
end
    
