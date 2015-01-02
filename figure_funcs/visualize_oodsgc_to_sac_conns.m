close all
c = 17080;
num_sacs = 10;
threshold = 150;

c_d = cell_data(c);
my_conts = double(c_d.contacts);
my_conts = my_conts(:,my_conts(2,:)>threshold);
is_sac = my_conts(1,:) > 70000 & my_conts(1,:) < 80000;
my_conts = my_conts(:,is_sac);

cont_sacs = unique(my_conts(1,:));
cont_counts = zeros(length(cont_sacs),1);

for n = 1:length(cont_sacs)
    cont_counts(n) = sum(my_conts(1,:) == cont_sacs(n));
end

[dummy sortord] = sort(cont_counts, 'descend');
    
cont_sacs = cont_sacs(sortord);

for n = 1:length(cont_sacs)
    
    cn = cont_sacs(n);
    c_d = cell_data(cn);
    is_me = my_conts(1,:) == cn;
    m = c_d.get_midpoint(true);
    figure; hold on
%     plot_cells(c, 1, .005, [1 1 1]*.7);
%     plot_cells(cn, 1, .005, [1 .6 .6]*.7);
    
    plot([m(2)*ones(1,sum(is_me)); my_conts(4,is_me)], [m(3)*ones(1,sum(is_me)); my_conts(5,is_me)]);

%     scatter([my_conts(4,is_me)', my_conts(5,is_me)', '*');
    
    if any(C.type.on_sac == cn);
        sac_type = 'ON';
    else
        sac_type = 'OFF';
    end
    
    title([num2str(cn) ' - ' sac_type]);
end



    