function sac2dsgc_vericose_contacts

    %right: 17080, 20213, 20210, 25005, 90001, 17080
    %left: 20179, 20210, 20239, 20245, 20254, 
    %up: 17161, 20233
    %down: 90002
    C = get_constants;
    
    colors = [1 0 0; 1 0 1; .8 .8 0; 0 0 0];
    dsgcs = [17080 20213 20210 25005 90001 17080 201789 20210 20239 20245 20542 17161 20233 90002];
    dsgc2color = [1 1 1 1 1 1 2 2 2 2 2 3 3 4];
    
    sacs = C.type.sure_off_sac;
    
    for n = 1:5;%length(sacs);
        c_d = cell_data(sacs(n));
        
        [conns, vericose_nodes] = detect_vericose_contacts(sacs(n), 1000, 4, 30000);
        if ~isempty(conns)

            soma_loc = c_d.get_midpoint;

            figure; hold all
            for k = 1:4
                is_my_color = dsgc2color == k;
                is_my_conn = false(1,size(conns,2));
                for l = find(is_my_color)
                    is_my_conn = is_my_conn | conns(1,:) == dsgcs(l);
                end
                sub_conn_loc = conns(4:5, is_my_conn);

                plot([soma_loc(2)*ones(1,size(sub_conn_loc,2)); sub_conn_loc(1,:)], ...
                    [soma_loc(3)*ones(1,size(sub_conn_loc,2)); sub_conn_loc(2,:)], ...
                    'lineWidth', 2, 'Color', colors(k,:));
            end
            
        end
    end

end