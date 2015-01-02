C = get_constants;

ts{1} = [C.type.t1 C.type.t2];
ts{2} = [C.type.t3a C.type.t3b C.type.t4];
ts{3} = C.type.sure_off_sac;
ts{6} = C.type.on_sac;
ts{4} = C.type.j;
ts{5} = C.type.minij;

strats = zeros(length(C.strat_x), 6);
for k = 1:5
    for c = ts{k}
        c_d = cell_data(c);
        s = c_d.stratification;
        strats(1:length(s),k) = strats(1:length(s),k) + s;
    end
    strats(:,k) = strats(:,k)/length(ts{k});
end

figure; hold all
for k = 1:5
    plot(C.strat_x(21:80), strats(21:80,k), 'lineWidth', 2);
end

legend('BC12', 'BC34', 'Off SAC', '"J"', 'mini-J')

