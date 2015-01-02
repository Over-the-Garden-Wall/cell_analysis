function [members, overlap] = cell_overlap(cell_nums)

    num_groups = length(cell_nums);

    num_cells = zeros(num_groups,1);
    for k = 1:num_groups
        num_cells(k) = length(cell_nums{k});
    end
    
    hulls = cell(num_groups,1);
    ind_hulls = cell(num_groups,1);
    
    ind_hulls{1} = cell(num_cells(1),1);
    for d = 1:num_groups
        for k = 1:num_cells(d)
            c_d = cell_data(cell_nums{d}(k));
            ind_hulls{d}{k} = cell(1,2);
            [ind_hulls{d}{k}{:}] = poly2cw(c_d.hull_2d(:,1), c_d.hull_2d(:,2));
        end
    end
    hulls{1} = ind_hulls{1};
    
    for d = 2:num_groups
        num_prior_hulls = length(hulls{d-1});
        
        hulls{d} = cell(num_prior_hulls*num_cells(d),1);

        for pd = 1:num_prior_hulls
            if ~isempty(hulls{d-1}{pd}) && ~isempty(hulls{d-1}{pd}{1})
                for cd = 1:num_cells(d)

                    hulls{d}{cd + (pd-1)*num_cells(d)} = cell(1,2);
                    [hulls{d}{cd + (pd-1)*num_cells(d)}{:}] = ...
                        polybool('intersection', ...
                        hulls{d-1}{pd}{1}, hulls{d-1}{pd}{2}, ...
                        ind_hulls{d}{cd}{1}, ind_hulls{d}{cd}{2});
                end
            
            end
        end
    end
    
    num_sets = prod(num_cells);
    orig_inds = cell(1,num_groups);
    [orig_inds{:}] = ind2sub(num_cells(end:-1:1)', 1:num_sets);
    members = zeros(num_sets,num_groups);
    for k = 1:num_groups
        members(:,k) = cell_nums{k}(orig_inds{end-k+1});
    end
    
    overlap = zeros(num_sets,1);
    for k = 1:num_sets
        if ~isempty(hulls{end}{k}) && ~isempty(hulls{end}{k}{1})
            overlap(k) = polyarea(hulls{end}{k}{1}, hulls{end}{k}{2});            
        end
    end
    
    [overlap, sortord] = sort(overlap, 'descend');
    members = members(sortord,:);
end
    
                
        