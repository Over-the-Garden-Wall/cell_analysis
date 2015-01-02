C = get_constants;
load(C.trans_loc);

for c = C.type.sure_off_sac
    disp(c);
    tic;
    cell_dat = cell_data(c);
    p = cell_dat.get_volume([1 1 1]'*[-Inf Inf]);
    [nodes edges root] = TEASAR(p, [16.5 16.5 25]*2, [250 250 250], 3, 3);
    for k = 1:size(nodes,1); 
        nodes(k,:) = apply_transform(T,nodes(k,:)); 
    end;
    save(['./skeletons/s' num2str(c) '.mat'], 'nodes', 'edges', 'root');
    toc
end