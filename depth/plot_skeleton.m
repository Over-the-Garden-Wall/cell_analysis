function plot_skeleton(skele, proj_plane, plot_color)

    if ~exist('plot_color','var') || isempty(plot_color);
        plot_color = [0 0 1];
    end
    
    
%     res = [16.5 16.5 25]*6;
    
%     res(proj_plane) = [];
    skele.nodes(:,proj_plane) = [];
%     
%     for d = 1:size(nodes,2)
%         nodes(:,d) = nodes(:,d)*res(d);
%     end
    
%     figure;
    
        
    for d = 1:size(skele.edges,1);
        n1 = skele.nodes(skele.edges(d,1),:);
        n2 = skele.nodes(skele.edges(d,2),:);
        
        plot([n1(1); n2(1)], [n1(2), n2(2)], 'Color', plot_color, 'LineWidth', 2);
        
        hold on
    end
end
        
    

    