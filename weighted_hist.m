function [hist_y, x] = weighted_hist(y, w, x)

    [dummy, x] = hist(y,x);
    
    hist_step = x(2)-x(1);
    h = hist_step/2;
    
    hist_y = zeros(size(dummy));
    for n = 1:length(x);
        hist_y(n) = sum(w(y>x(n)-h & y <= x(n)+h));
    end
    
    plot(x, hist_y, 'lineWidth', 2);
end