function sac_to_sac_connections
tic
C = get_constants;

info_struct = contact_info(C.type.sure_off_sac, C.type.sure_off_sac);

cell_midsa = info_struct.soma_loc(info_struct.cell_ids(:,1),2:3);
cell_midsb = info_struct.soma_loc(info_struct.cell_ids(:,2),2:3);

dista = sqrt((cell_midsa(:,1) - info_struct.contact_loc(:,2)).^2 + ...
    (cell_midsa(:,2) - info_struct.contact_loc(:,3)).^2);
distb = sqrt((cell_midsb(:,1) - info_struct.contact_loc(:,2)).^2 + ...
    (cell_midsb(:,2) - info_struct.contact_loc(:,3)).^2);

angles = acos(-((cell_midsa(:,1)-info_struct.contact_loc(:,2)) .* ...
    (cell_midsb(:,1)-info_struct.contact_loc(:,2)) + ...
    (cell_midsa(:,2)-info_struct.contact_loc(:,3)) .* ...
    (cell_midsb(:,2)-info_struct.contact_loc(:,3))) ./ dista ./ distb);

angles = ceil(angles/2/pi*360);
cont_area = zeros(180,1);
for x = 1:180
    cont_area(x) = mean(info_struct.contact_size(angles==x));
end

figure; plot(cont_area)

% hist_denom = hist_denom+hist_denom(end:-1:1);
% hist_denom = hist_denom(1:length(hist_data));
% figure; plot(hist_data./hist_denom');



