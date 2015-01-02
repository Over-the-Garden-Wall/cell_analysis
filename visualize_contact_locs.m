C = get_constants;

sac_num = C.type.sure_off_sac(1);
cont_thresh = 1;


figure; hold all;

plot_cells(sac_num, 1, .01, [0 0 0]);

cell_dat = cell_data(sac_num);

conts = cell_dat.contacts;
conts(:,conts(2,:)<cont_thresh) = [];
num_conts = size(conts,2);

is_type_2 = false(num_conts,1);
is_type_3a = false(num_conts,1);

for n = 1:num_conts;
    is_type_2(n) = any(conts(1,n)==C.type.t2);
    is_type_3a(n) = any(conts(1,n)==C.type.t3a);
end

scatter(conts(4,is_type_2),conts(5,is_type_2),'*');
scatter(conts(4,is_type_3a),conts(5,is_type_3a),'*');

    