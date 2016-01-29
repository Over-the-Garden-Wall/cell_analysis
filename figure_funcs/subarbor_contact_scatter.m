%analyze sac subarbors

C = get_constants;
cns = C.type.on_sac_subarbor;
num_cells = length(cns);

contacting_types = {'BC7', 'BC5t'};
num_types = length(contacting_types);

contact_means = zeros(num_cells, num_types);

max_cell_num = 99999;
ID2type = zeros(max_cell_num,1);
for t = 1:num_types
    ID2type(C.type.(contacting_types{t})) = t;
end
    

for n = 1:num_cells
    c_d = cell_data(cns(n));
    
    my_conts = c_d.contacts;
    my_soma = c_d.get_midpoint(true);
    
    for t = 1:num_types
        sub_conts = my_conts(2:end, ID2type(my_conts(1,:))==t);
        
        cont_d = sqrt((my_soma(2)-sub_conts(3,:)).^2 + (my_soma(3)-sub_conts(4,:)).^2);
        cont_w = sub_conts(1,:) / sum(sub_conts(1,:));
        
        contact_means(n,t) = sum(cont_d.*cont_w);
    end
end

figure; hold all

for t = 1:num_types
    scatter(1:num_cells, contact_means(:,t)', '*');
end

figure; hold all
my_clr = [.7 .7 0];
old_cn = 0;
for n = 1:num_cells
    my_cn = cns(n);
    if my_cn ~= old_cn + 100000
        my_clr = 1-my_clr;
    end
    old_cn = my_cn;
    scatter(n, contact_means(n,2)-contact_means(n,1), '*', 'markerEdgeColor', my_clr);
end


sac_nums = [];
sac_means = [];
for n = 1:num_cells
    my_sac_num = mod(cns(n), 100000);
    my_k = (cns(n) - my_sac_num) / 100000;
   
    if my_k == 1
        sac_nums(end+1) = my_sac_num;
        sac_means(end+1) = contact_means(n,2)-contact_means(n,1);
    else
        sac_means(end) = (contact_means(n,2)-contact_means(n,1))/k + sac_means(end)/k*(k-1);
    end
end

[a, b] = sort(sac_means);

figure; hold all
for n = 1:length(sac_nums);
    c_d = cell_data(sac_nums(n));
    sac_mid = c_d.get_midpoint(true);
    
    text(sac_mid(2), sac_mid(3), num2str(b(n)));
end
    
