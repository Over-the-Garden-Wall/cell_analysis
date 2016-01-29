C = get_constants;
load(C.trans_loc);

% types = {'on_sac', 'sure_off_sac'};
% types = {'t1','t2','t3a','t3b','t4','t5w','t5h','t5l','t6','t7','t89', 'xbc','tRBC'};
% types = {'biGC1'};
types = {'gc7i', 'gc37'};
dsmp = [2 2 2];
%stopped on 17097
for t = 1:length(types)
    cns = C.type.(types{t});
%     cns = cns
    if t == 1
        warning('hackzor');
        a = find(cns==17097);
        cns(1:a-1) = [];
    end
    for c = cns
        
            
%         for c = 90001
    disp(c);
    tic;
    cell_dat = cell_data(c);
    
    p = cell_dat.get_volume([1 1 1]'*[-Inf Inf]);
    [nodes edges root node_diameter] = TEASAR(p, dsmp, 3, 3);
    node_diameter = node_diameter * mean(C.res .* dsmp);
    save([C.skele_untrans_dir '/s' num2str(c) '.mat'], 'nodes', 'edges', 'root', 'node_diameter');
    for k = 1:size(nodes,1); 
        nodes(k,:) = apply_transform(T,nodes(k,:)); 
    end;
    save([C.skele_dir '/s' num2str(c) '.mat'], 'nodes', 'edges', 'root', 'node_diameter');
    toc
end

end