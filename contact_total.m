function cont = contact_total(cell_ids1, cell_ids2)
    if ~iscell(cell_ids1)
        cell_ids1 = {cell_ids1};
    end
    if ~iscell(cell_ids2)
        cell_ids2 = {cell_ids2};
    end
    
    cont = zeros(length(cell_ids1), length(cell_ids2));
    
    for k = 1:length(cell_ids1)        
        for c = cell_ids1{k}
            c_d = cell_data(c);
            conts = double(c_d.contacts);
            
            for kk = 1:length(cell_ids2)
                for c2 = cell_ids2{kk}
%                     if c_d.contact_map.isKey(c2)
%                         cont(k,kk) = cont(k,kk) + c_d.contact_map(c2);
%                     end
                    cont(k,kk) = cont(k,kk) + sum(conts(2,conts(1,:)==c2));
                end
            end
        end
    end
    
end
                    