C = get_constants;

cell_nums = C.type.t5;

num_cells = length(cell_nums);
G = zeros(num_cells);

for mk = 1:num_cells;
    for nk = mk+1:num_cells;
        c_d1 = cell_data(cell_nums(mk));
        c_d2 = cell_data(cell_nums(nk));
        
        h1 = [];
        [h1(:,1), h1(:,2)] = poly2cw(c_d1.hull_2d(:,1), c_d1.hull_2d(:,2));
        
        h2 = [];
        [h2(:,1), h2(:,2)] = poly2cw(c_d2.hull_2d(:,1), c_d2.hull_2d(:,2));
        
        int_hull = [];
        try
            [int_hull(:,1), int_hull(:,2)] = polybool('intersection',h1(:,1), h1(:,2), h2(:,1), h2(:,2));
            G(mk,nk) = polyarea(int_hull(:,1), int_hull(:,2));
        catch ME
            
        end
    end
end
   
G = G+G';

G = G/max(G(:));

misfits = sum(G)==0;

misfit_nums = cell_nums(misfits);
cell_nums(misfits) = [];
G = G(~misfits,:);
G = G(:,~misfits);

num_nonzeros = sum(G(:)~=0);
fid = fopen('./biqmac_out.txt','w');
fprintf(fid, '%d %d\n', length(cell_nums), num_nonzeros);
for y = 1:size(G,1)
    for x = 1:size(G,2)
        if G(y,x)~=0
            fprintf(fid, '%d %d %d\n', y, x, round(G(y,x)*10000));
        end
    end
end
fclose(fid);

% 
% Warning: multiple edges! Edge weights have been summed up.
% The input graph has 71 nodes and 386 edges, edge weights are integers.
% 
% Total number of branch-and-bound nodes: 107
% Total time: 373.31 sec
% cut value: 1029412.000000
% one side of the cut:
%  2 3 5 6 7 9 10 12 13 14 18 19 20 21 24 25 27 28 32 34 37 38 41 44 47 50 52 54 56 58 60 61 63 66 68

% Warning: multiple edges! Edge weights have been summed up.
% The input graph has 62 nodes and 322 edges, edge weights are integers.
% 
% Total number of branch-and-bound nodes: 49
% Total time: 111.30 sec
% cut value: 941814.000000
% one side of the cut:
%  2 4 6 7 11 12 15 16 22 23 26 28 30 33 35 36 37 39 40 41 42 44 47 49 50 51 56 59 61 62