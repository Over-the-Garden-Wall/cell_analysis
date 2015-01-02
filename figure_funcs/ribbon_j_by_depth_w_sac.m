percs = [1/3 2/3];
body_cutoff = 85;



C = get_constants;

num_cells = length(C.type.j);

num_percs = length(percs)+1;

strats = cell(num_percs,1);
mean_strat = cell(num_percs,1);
strat_ste = cell(num_percs,1);
strat_x = cell(num_percs,1);

min_length = 200;
for n = 1:3;
    strats{n} = zeros(min_length,num_cells);
    mean_strat{n} = zeros(min_length,1);
    strat_ste{n} = zeros(min_length,1);
end



for n = 1:num_cells;

    cell_dat = cell_data(C.type.j(n));
    p = cell_dat.get_point_count;
    
    p(body_cutoff:end,:,:) = [];
    
    distance_distribution = sum(sum(double(p),1),3);
    cum_dist = cumsum(distance_distribution);
    
    min_length = min(min_length, size(p,1));
    
    
    for k = 1:3
        if k == 1
            lower_bound = 0;
        else
            lower_bound = percs(k-1);
        end
        if k == num_percs
            upper_bound = 1;
        else
            upper_bound = percs(k);
        end
                
        lower_range = find(cum_dist > cum_dist(end)*lower_bound,1,'first');
        upper_range = find(cum_dist <= cum_dist(end)*upper_bound,1,'last');   

        strats{k}(1:min_length,n) = sum(sum(p(1:min_length,lower_range:upper_range,:),3),2);
        strats{k}(:,n) = strats{k}(:,n)/sum(strats{k}(:,n));
    end
end

for k = 1:num_percs
    strats{k} = strats{k}(1:min_length,:);
    strat_ste{k} = std(strats{k},[],2)/sqrt(num_cells-1);
    strat_ste{k}(isnan(strat_ste{k})) = 0;
    mean_strat{k} = mean(strats{k},2);
    
    strat_x{k} = C.strat_x(1:min_length)';
end
    
    


[a b] = ribbon_plotn(strat_x, mean_strat, strat_ste);


% t12 = [C.type.t1 C.type.t2];
% t34 = [C.type.t3a C.type.t3b C.type.t4];
% 
% t_strats = zeros(200,2);
% for n = t12
%     cell_dat = cell_data(n);
%     cs = cell_dat.stratification;
%     t_strats(1:length(cs),1) = t_strats(1:length(cs),1) + cs;
% end
% for n = t34
%     cell_dat = cell_data(n);
%     cs = cell_dat.stratification;
%     t_strats(1:length(cs),2) = t_strats(1:length(cs),2) + cs;
% end
% 
% t_strats = t_strats(1:min_length,:);
% t_strats(:,1) = t_strats(:,1)/sum(t_strats(:,1));
% t_strats(:,2) = t_strats(:,2)/sum(t_strats(:,2));
% 
% 
% l = plot(x, t_strats(:,1)', x, t_strats(:,2)');
% legend([b; l], {'proximal','distal', 'type 1/2', 'type 3/4'});

off_sac = C.type.off_sac;

t_strats = zeros(200,1);
for n = off_sac
    cell_dat = cell_data(n);
    cs = cell_dat.stratification;
    t_strats(1:length(cs),1) = t_strats(1:length(cs),1) + cs;
end

t_strats = t_strats(1:min_length,:);
t_strats(:,1) = t_strats(:,1)/sum(t_strats(:,1));


l = plot(strat_x{k}, t_strats(:,1)');
legend([b; l], {'proximal', 'mid', 'distal', 'off_sac'});


