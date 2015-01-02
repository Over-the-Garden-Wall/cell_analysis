C = get_constants;

bip_types = {'off_sac'};


    [mean_loc total_contact] = get_mean_axial_contacts(10010, C.type.(bip_types{k}), true);
    [mean_loc, sort_inds] = sort(mean_loc{k});
    total_contact{k} = total_contact{k}(sort_inds);
% 
% for k = 1:5;
%     d{k} = sum(sum(d{k},1),3);
%     dc{k} =sum(sum(dc{k},1),3);    
% end
% 
% stepper = 10;
% for k = 1:5;
%     for bin = 1:ceil(length(d{k})/stepper)
%         b{k}(bin) = sum(d{k}(((bin-1)*stepper + 1): min(bin*stepper,end)));
%         bc{k}(bin) = sum(dc{k}(((bin-1)*stepper + 1): min(bin*stepper,end)));
%     end
%     
%     b{k} = b{k}./bc{k};
% %     b{k} = b{k} / length(C.type.(bip_types{k}));
%     b{k}(isnan(d{k})) = 0;
% end




figure; hold all
% for k = 1:5;
%     plot((C.axial_x_min:stepper:length(d{k}) + C.axial_x_min - 1) + stepper/2, b{k});
    scatter(mean_loc, total_contact, '*');

    
% end







% figure; hold all
% for k = 1:5;
%     plot((C.axial_x_min:stepper:length(d{k}) + C.axial_x_min - 1) + stepper/2, bc{k});
% 
%     
% end

