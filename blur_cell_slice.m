function im = blur_cell_slice(cell_num, blur_window_size, slice_min, varargin)

    p = inputParser;
    
    p.addRequired('cell_num', @(x)isnumeric(x) && length(x)==1);
    p.addRequired('blur_window_size', @(x)isnumeric(x) && length(x)==1);
    p.addRequired('slice_min', @(x)isnumeric(x) && length(x)==1);
    
    p.addOptional('slice_max', slice_min+1, @(x)isnumeric(x) && length(x)==1);
    p.addOptional('im_max', [], @(x)isnumeric(x) && length(x)==2);
    
    p.parse(cell_num, blur_window_size, slice_min, varargin{:});
    
    s = p.Results;

    C = get_constants;
    cx = cell_data(cell_num);
    
    px = cx.get_surface;
    px(:,1) = C.f(px(:,1));
    px = px(px(:,1)>=slice_min & px(:,1)<s.slice_max,2:3);
    
    px = ceil(px/1000);
    
    
    if isempty(s.im_max)
        im = zeros(max(px(:,1)), max(px(:,2)));
    else
        im = zeros(s.im_max);
    end
    
    for k = 1:size(px,1);
        im(px(k,1),px(k,2)) = im(px(k,1),px(k,2))+1;
    end
    
%     K = gausswin(blur_window_size);
%     K = K*K';
%     K = K/sum(K(:));
    K = ones(blur_window_size);
    K = K/sum(K(:));

    
    im = conv2(im,K,'same');
    
end
    