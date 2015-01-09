C = get_constants;

fns = dir(C.raw_point_dir);

cell_nums = zeros(length(fns),1);
for n = 1:length(fns)
    
    unders = strfind(fns(n).name,'_');
    if length(unders)<2
        cell_nums(n) = 0;
    else
        cell_nums(n) = str2double(fns(n).name(unders(1)+1:unders(2)-1));        
    end
    
end



cell_nums = unique(cell_nums);
cell_nums(cell_nums==0) = [];

for n = 1:length(cell_nums)
    fn = [C.surface_point_dir 'cell_' num2str(cell_nums(n)) '_surface.mat'];
    
    if ~exist(fn,'file')
        disp(fn)
        surface_points = find_surface_points(cell_nums(n));
    
        save(fn,'surface_points');
    end
end