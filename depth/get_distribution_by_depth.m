function [x y y2] = get_distribution_by_depth(target_cell, contacting_cells, R, Q)

res = [16.5 16.5 25];


new_R = double(R);
for k=1:3
    new_R(3+k,:) = new_R(3+k,:)*res(k);
end
new_R(4:6,:) = Q*new_R(4:6,:);


is_good = false(1,size(R,2));    
for k = contacting_cells; 
        is_good = is_good | (R(1,:)==target_cell & R(2,:)==k);
        is_good = is_good | (R(2,:)==target_cell & R(1,:)==k);
        
end

new_R = new_R([3 4],is_good)';

[dummy, sort_ind] = sort(new_R(:,2));
new_R = new_R(sort_ind,:);

x_min = min(new_R(:,2));
x_max = max(new_R(:,2));

x = floor(x_min):ceil(x_max);
new_R(:,2) = round(new_R(:,2));

y = zeros(1,length(x));

y(new_R(:,2)-floor(x_min)+1) = new_R(:,1);

y2 = conv(y, gausswin(1000)/sum(gausswin(1000)), 'same');

% x = x -;

f = inline('-.006*(x-2.7*10^4)+30');
x = f(x);


figure;subplot(1,2,1); plot(x,y2)
subplot(1,2,2); plot(x,y)

end