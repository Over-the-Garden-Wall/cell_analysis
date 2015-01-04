function setup_path
    root_path = pwd;

    addpath(root_path);
    addpath([root_path '/scripts/']);
    addpath([root_path '/depth/']);
    addpath([root_path '/plane_fitting/']);
    addpath([root_path '/cell_data_funcs/']);

    addpath([root_path './tree_analysis/']);
%     
%     addpath('~greenem/stratification/matlab_bgl/');
%     addpath('~greenem/stratification/DSDP5.8_64/matlab/');
end
    