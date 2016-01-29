function rasterize_figure(ax_h)

    axis_properties = get(ax_h);
    
    fig_h = axis_properties.Parent;
    
    fig_c = get(fig_h, 'Color');
    set(fig_h, 'Color', [1 1 1]);
    
    
    set(ax_h, 'Visible', 'off');
    F = getframe(ax_h);
    
    set(fig_h, 'Color', fig_c);
    
    M = F.cdata;
    cla(ax_h);
    axes(ax_h);
    im_h = image(M);
    set(ax_h, 'Visible', 'on');    
    
%     im_h = imagesc(M(end:-1:1,:,:), [0 255]);
    uistack(im_h, 'bottom')
    new_ax = gca;
    new_xtick = (axis_properties.XTick - axis_properties.XLim(1)) / (axis_properties.XLim(2) - axis_properties.XLim(1)) * (size(M,2)-1) + 1;
    new_ytick = (axis_properties.YTick - axis_properties.YLim(1)) / (axis_properties.YLim(2) - axis_properties.YLim(1)) * (size(M,1)-1) + 1;
    
    
    
    set(new_ax, 'XTick', new_xtick, 'XTickLabel', axis_properties.XTickLabel, 'XLim', [1 size(M,2)],...
        'YTick', new_ytick, 'YTickLabel', axis_properties.YTickLabel, 'YLim', [1 size(M,1)]);
    
%     if strcmp(axis_properties.YDir, 'normal')
%         set(new_ax, 'YDir', 'reverse');
%     else
%         set(new_ax, 'YDir', 'normal');
%     end
    
%     close(fig_h)
end
        