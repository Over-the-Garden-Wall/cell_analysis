
C = get_constants;

types = {'xbc', 't5w', 't5', 't6', 't7', 't89', 'rbc'};

ranges = {[56 58], [49.5 52.5; 59 62], [52 61], [56 60; 69 73], [65 73], [75 85], [85 95]};


all_cells = C.type.on_bc;
perc_vals = zeros(length(all_cells),length(types));

for k = 1:length(types)
   
    perc_vals(:,k) = percent_strat_within_range(all_cells, ranges{k}, [40 100]);
end