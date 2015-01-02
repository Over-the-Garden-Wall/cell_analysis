batch_find_surface;

load('./conns.mat');
load('./T.mat');

C = get_constants;

transform_points('./surface_points',C.point_dir, T, true);
conns = transform_conns(conns,T);
save(C.conn_loc,'conns');