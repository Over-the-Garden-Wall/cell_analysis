C = get_constants;

CELL_ID_MAX = 100000;

type_names = {'BC5i', 'BC5t', 'BC5o'};
BC5s = [];
BC5_labels = [];

for t = 1:length(type_names);
    BC5s = [BC5s, C.type.(type_names{t})];
    BC5_labels(end+1:length(BC5s)) = t;
end

num_BCs = length(BC5s);


GCs = C.type.gc63;
num_GCs = length(GCs);

GC_ID = zeros(CELL_ID_MAX, 1);
GC_ID(GCs) = 1:length(GCs);


BC_SA = zeros(num_BCs,1);
BC_GC_conn = zeros(num_BCs, num_GCs);

for bc = 1:length(BC5s);
    c_d = cell_data(BC5s(bc));
    db_conts = double(c_d.contacts);
    
    BC_SA(bc) = c_d.SA;
    
    for gc = 1:num_GCs
        BC_GC_conn(bc, gc) = sum(db_conts(2, GC_ID(db_conts(1,:))==gc));
    end
end
    
    