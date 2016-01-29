function supplementary_table(types)

    C = get_constants;
    num_types = length(types);
    
    hull_areas = cell(num_types,1);
    union_area = zeros(num_types,1);
        
    
    for t = 1:num_types
        cns = C.type.(types{t});
        total_hull = cell(1,2);
        h = cell(1,2);
        
        hull_areas{t} = zeros(length(cns),1);
        
        for cn = 1:length(cns);
            
            c_d = cell_data(cns(cn));
            [h{:}] = poly2cw(c_d.hull_2d(:,1), c_d.hull_2d(:,2));
            [total_hull{:}] = polybool('union', h{1}, h{2}, total_hull{1}, total_hull{2});
            
            hull_areas{t}(cn) = c_d.hull_area;
        end
        union_area(t) = poly_area([total_hull{1}, total_hull{2}]);
    end
    
        fprintf('\n');
    for t = 1:num_types
        out_line = {types{t}, ...
            num2str(length(C.type.(types{t}))), ...
            [num2str(round(mean(hull_areas{t}/10^6))) '+-' num2str(round(std(hull_areas{t}/10^6)))], ...
            num2str(round(length(C.type.(types{t})) / union_area(t) * 10^12)), ...
            num2str(round(sum(hull_areas{t}) / union_area(t) * 100 - 100))};
        
        for n = 1:length(out_line)
            fprintf('%s\t', out_line{n});
        end
        fprintf('\n');
    end
end