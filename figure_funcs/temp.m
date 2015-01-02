k = dir('~/stratification/point_data');


for n = 1:length(fns)
    if k(n).name(1)=='o'
        disp(k(n).name)
        p_num = ['10' k(n).name(11:15)];
        movefile(['~/stratification/point_data/' k(n).name], ['~/stratification/point_data/points_' p_num k(n).name(16:end)]);
    end
end