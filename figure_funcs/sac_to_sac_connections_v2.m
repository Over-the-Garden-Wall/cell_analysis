function sac_to_sac_connections
tic
C = get_constants;

info_struct = contact_info(C.type.on_sac, C.type.on_sac);

cell_mids = cell(2,1);
abs_vec = cell(2,1);
dist = cell(2,1);

for k = 1:2
    cell_mids{k} = info_struct.soma_loc(info_struct.cell_ids(:,k),2:3);
    abs_vec{k} = cell_mids{k} - info_struct.contact_loc(:,2:3);
    dist{k} = sqrt(sum(abs_vec{k}.^2,2));
end

angles = acos(sum(abs_vec{1}.*abs_vec{2}, 2) ./ dist{1} ./ dist{2});





x_bounds = [0 3.9*10^5];
y_bounds = [10^4 3.4*10^5];


rads = (1:360) / 180 * pi;
ucircle = [cos(rads') sin(rads')];

deg_count = zeros(360);

for k = 1:size(angles,1)
    if mod(k, 10000) == 1
        disp(['working on contact #' num2str(k) ' of ' num2str(size(angles,1))]);
        toc
        tic
    end
    
    circax = ucircle(:,1)*dista(k) + cell_midsa(k,1);
    circay = ucircle(:,2)*dista(k) + cell_midsa(k,2);
    
    circbx = ucircle(:,1)*distb(k) + cell_midsb(k,1);
    circby = ucircle(:,2)*distb(k) + cell_midsb(k,2);
    
    is_ina = double(circax > x_bounds(1) & circax < x_bounds(2) & circay > y_bounds(1) & circay < y_bounds(2));
    is_inb = double(circbx > x_bounds(1) & circbx < x_bounds(2) & circby > y_bounds(1) & circby < y_bounds(2));
    
    deg_count = deg_count + is_ina*is_inb';
end

hist_data = hist(angles, rads(1:length(rads)/2)-(rads(2)-rads(1))/2);
hist_denom = zeros(360,1);

for x = 1:360
    for y = 1:360
        theta = max(y,x) - min(y,x) + 1;
        hist_denom(theta) = hist_denom(theta) + deg_count(x,y);
    end
end

hist_denom = hist_denom+hist_denom(end:-1:1);
hist_denom = hist_denom(1:length(hist_data));
figure; plot(hist_data./hist_denom');



