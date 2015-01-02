function histo_data = weighted_histogram(data, weights, bin_spec)

    if ~exist('bin_spec','var') || isempty(bin_spec)
        bin_spec = 10;
    end
    
    if length(bin_spec) == 1
        bin_step = (max(data) - min(data)) / bin_spec;
        bin_spec = min(data) + (.5:bin_spec)*bin_step;
    end
    
    num_bins = length(bin_spec);
    
    bins = [-Inf sort(bin_spec) Inf];    
    
    histo_data = zeros(num_bins,1);
    for n = 1:num_bins
        is_this_bin = data >= (bins(n) + bins(n+1))/2 & data < (bins(n+1) + bins(n+2))/2;
        histo_data(n) = sum(weights(is_this_bin));
    end
    
    plot(histo_data);
end