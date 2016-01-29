batch_find_surface;

C = get_constants;

load(C.raw_conn_loc);
load(C.trans_loc);

transform_points(C.surface_point_dir,C.point_dir, T, true);
conns = transform_conns(conns,T);
save(C.conn_loc,'conns');

make_stuff_for_jinseop;