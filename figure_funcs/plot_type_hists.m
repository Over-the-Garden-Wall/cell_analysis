
% plot_point = 25;
types = {'t1','t2'};

hist_data = cell(length(types),1);
for k = 1:length(types)
    hist_data{k} = perc_dat{k}(:,75) - perc_dat{k}(:,25);
end

maxval = -Inf; minval = Inf; 

for k = 1:length(types); 
    
    minval = min([hist_data{k}; minval]); 
    maxval = max([hist_data{k};maxval]); 
end


r = floor(minval):ceil(maxval);

hists = zeros(length(r),length(types));
for k = 1:length(types); 
    hists(:,k) = hist(hist_data{k},r); 
end

h = area(r, hists);
for k = 1:length(types); 
    set(h(k), 'FaceColor', C.colormap(k,:))
end
prep_figure(gcf, gca, 'xlabel', ['75^t^h - 25^t^h percentile'], 'legend', {'BC1','BC2','BC3a','BC3b','BC4'});