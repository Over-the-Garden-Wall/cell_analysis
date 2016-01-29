function save_fig(fn, fig_h)

    if ~exist('fig_h','var') || isempty(fig_h)
        fig_h = gcf;
    end
    set(fig_h, 'PaperPositionMode', 'auto');
    print(fn, '-depsc', '-r0');
    
end