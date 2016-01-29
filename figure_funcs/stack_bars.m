function stack_bars(my_axis, clrs)

    ax_children = get(my_axis, 'Children');
       
    clr_counter = 0;
    for n = 1:length(ax_children)
        all_props = get(ax_children(n));
        if isfield(all_props, 'Vertices');
            clr_counter = clr_counter+1;
            my_verts = all_props.Vertices;

            if ~exist('bar_bottoms', 'var')
                bar_bottoms = zeros(floor(size(my_verts,1)/5),1);
                bar_mids = my_verts(3:5:end,1);
                bar_width = my_verts(6,1) - bar_mids(1);

            end

            bar_heights = my_verts(3:5:end,2);

            my_verts((1:length(bar_mids))*5 - 4,1) = (bar_mids - bar_width/2);
            my_verts((1:length(bar_mids))*5 - 3,1) = (bar_mids - bar_width/2);
            my_verts((1:length(bar_mids))*5 - 2,1) = (bar_mids - bar_width/2);
            my_verts((1:length(bar_mids))*5 - 1,1) = (bar_mids + bar_width/2);
            my_verts((1:length(bar_mids))*5,1) = (bar_mids + bar_width/2);

            my_verts((1:length(bar_mids))*5 - 4,2) = bar_bottoms;
            my_verts((1:length(bar_mids))*5 - 3,2) = bar_bottoms;
            my_verts((1:length(bar_mids))*5 - 2,2) = bar_bottoms + bar_heights;
            my_verts((1:length(bar_mids))*5 - 1,2) = bar_bottoms + bar_heights;
            my_verts((1:length(bar_mids))*5,2) = bar_bottoms;

            my_verts(length(bar_mids)*5+1:end,1) = my_verts(length(bar_mids)*5,1);
            my_verts(length(bar_mids)*5+1:end,2) = my_verts(length(bar_mids)*5,2);

            bar_bottoms = bar_heights + bar_bottoms;

    %         set(ax_children(n), 'Vertices', my_verts)
            if exist('clrs', 'var')
                set(ax_children(n), 'faceColor', clrs(clr_counter,:));            
            end
            set(ax_children(n), 'Vertices', my_verts, 'edgeColor', get(ax_children(n), 'faceColor'));
        
        end
    end
    
    
end
    
    