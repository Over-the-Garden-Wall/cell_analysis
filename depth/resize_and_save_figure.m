function resize_and_save_figure(h, save_name)

    set(h,'Position',[1 1 1080 800]);
    saveas(h,save_name,'tif')

end