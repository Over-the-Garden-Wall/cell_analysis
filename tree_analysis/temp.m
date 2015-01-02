C = get_constants;
types = {'sure_off_sac', 't2','t3a'};

figure; hold all;
strats = zeros(200,3); 
for k = 1:3; 
    for c = C.type.(types{k}); 
        c_d = cell_data(c); 
        s = c_d.stratification; 
        strats(1:length(s),k) = strats(1:length(s),k) + s; 
    end; 
    strats(:,k) = strats(:,k) / length(C.type.(types{k})); 
    plot_lim = find(strats(:,k)~=0,1,'last');
    plot(C.strat_x(21:plot_lim), strats(21:plot_lim,k), 'LineWidth', 2);
end

legend(types);


figure; hold all;

for k = 2:3;
    plot_vals = strats(:,1).*strats(:,k);
    plot_lim = find(plot_vals~=0,1,'last');
    plot(C.strat_x(21:plot_lim), plot_vals(21:plot_lim), 'LineWidth', 2);
end

legend(types(2:3));



figure; hold all

cont_strats = zeros(200,2);
for k = 1:2
    for c = C.type.sure_off_sac
        c_d = cell_data(c);
        conts = c_d.contacts;
        for t = 1:size(conts,2)
            if any(conts(1,t) == C.type.(types{k+1}))
                c_loc = round(C.f(conts(3,t)))-C.strat_x(1)+1;
                cont_strats(c_loc,k) = cont_strats(c_loc,k) + conts(2,t);
            end
        end
    end
    
    plot_lim = find(cont_strats(:,k)~=0,1,'last');
    plot(C.strat_x(21:plot_lim), cont_strats(21:plot_lim,k), 'LineWidth', 2);
end

legend(types(2:3));

