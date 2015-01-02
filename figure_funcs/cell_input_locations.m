function type_data = cell_input_locations(ref_cell, types, use_axis)

    if ~exist('use_axis','var') || isempty(use_axis)
        use_axis = false;
    end

    DATA_MAX = 10000;
    C = get_constants;
    
    
    num_types = length(types);
    c_d = cell_data(ref_cell);
    
    type_data = cell(num_types,1);
    
    soma_loc = c_d.get_midpoint(true);
    
    
    for k = 1:num_types
        type_data{k} = zeros(DATA_MAX, 5);
        dcount = 1;
        for n = 1:size(c_d.contacts,2);
            if any(c_d.contacts(1,n)==C.type.(types{k}))
                type_data{k}(dcount,1:4) = double(c_d.contacts(2:5,n));
                dcount = dcount+1;
            end
        end
        
        if use_axis
            d = (type_data{k}(:,3)-soma_loc(2))*c_d.dist_axis(1) + (type_data{k}(:,4)-soma_loc(3))*c_d.dist_axis(2);
        else
            d = sqrt((type_data{k}(:,3)-soma_loc(2)).^2 + (type_data{k}(:,4)-soma_loc(3)).^2);
        end
        type_data{k}(:,5) = d;
        type_data{k}(dcount:end,:) = [];
    end
    
%     figure; hold all
%     for k = 1:num_types
%         scatter(type_data{k}(:,3), type_data{k}(:,4));
%     end
%     legend(types);
%     
%     
%     figure; hold all
%     for k = 1:num_types
%         
%         scatter(type_data{k}(:,5),type_data{k}(:,1));
%     end
%     legend(types);
%     
%     
%     figure; hold all
%     for k = 1:num_types
%         
%         
%         K = gausswin(10000);
%         K = K/sum(K);
%         
%         x = round(min(type_data{k}(:,5))):round(max(type_data{k}(:,5)));
%         y = zeros(size(x));
%         for n = 1:length(type_data{k}(:,5));
%             d_loc = round(type_data{k}(n,5)) - round(min(type_data{k}(:,5))) + 1;
% %             y(d_loc) = y(d_loc) + type_data{k}(n,1);
%             y(d_loc) = y(d_loc) + 1;
%         end
%         
%         y = conv(y,K,'same');
%         plot(x,y);
%     end
%     
%     legend(types);
end
    