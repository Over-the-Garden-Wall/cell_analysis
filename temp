
y = zeros(130,5);
for k = 1:5;
    for x = 1:130;
        is_valid = contact_list{k}(:,1) == x;
        if any(is_valid)
            y(x,k) = mean(contact_list{k}(is_valid,2));
        end
    end
end