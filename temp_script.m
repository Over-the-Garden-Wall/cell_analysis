
%sac skeleton-fixing script
C = get_constants;
load(C.trans_loc);

fns = dir(C.skele_untrans_dir);

for n = 1:length(fns)
    fn = fns(n).name;
    if fn(1) == 's'
        load([C.skele_untrans_dir, '/', fn]);
        nodes = apply_transform(T, nodes);
        save([C.skele_dir, '/', fn], 'nodes', 'edges', 'root', 'node_diameter');
    end
end