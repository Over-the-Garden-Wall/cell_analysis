function out_sig = dsgc_model(I, v)

    gc_d = load('~/data/stratification/model_dump.mat');
    dsgc_radius = 100;
    
    a = find(sum(gc_d.dist_hist(end:-1:1,:),2), 1, 'first');
    
    neuron_irf(:,1) = gradient(lognpdf(0:.001:1, .1, 1.25));
    in_sign(1) = 1;
    neuron_irf(:,2) = gradient(lognpdf(0:.001:1, .1, 1.25));
    in_sign(2) = -1;
    
    
    sac_max_length = size(gc_d.dist_hist, 1) - a + 1;
    
    
    
    im_sz = dsgc_radius*2 + sac_max_length;
    if mod(im_sz, 2) == 0
        im_sz = im_sz + 1;
    end
    im_center = (im_sz+1) / 2;
    
    
    rf_im = zeros(im_sz, im_sz, 2, 2);
    
    [X, Y] = meshgrid( (1:im_sz)-im_center, (1:im_sz)-im_center );
    R = sqrt(X.^2 + Y.^2);
    A = atan2(Y, X) + pi;
    
    for l = 1:2
        sac_dend = gc_d.dist_hist(:,l) / sum(gc_d.dist_hist(:,l));
        sac_dend = cumsum(sac_dend(end:-1:1));
        sac_dend = sac_dend(end:-1:1);
        
        rf_im(:,:,1, l) = R < dsgc_radius;
    
        sac_im = zeros(im_sz, im_sz);
        for y = 1:im_sz
            for x = 1:im_sz                
                r = round(R(y,x));
                a = ceil(A(y,x) / 2 / pi * size(gc_d.angle_distribution, 1));
                if r > 0 && r < size(gc_d.dist_hist, 1)
                    sac_im(y, x) = sac_dend(r) * gc_d.angle_distribution(a, l);
                end
            end
        end
        rf_im(:,:,2, l) = conv2(rf_im(:,:,1, l), sac_im, 'same');
                        
    end
    
    for l = 1:2
        for k = 1:size(rf_im, 3)
            rf_im(:,:,k,l) = rf_im(:,:,k,l) / abs(sum(sum(rf_im(:,:,k,l))));
        end
    end
    
    sz_diff = size(I,1) - size(rf_im, 1);
    if sz_diff > 0
        I = I(1 + floor(sz_diff/2):end - ceil(sz_diff/2),:);
    elseif sz_diff < 0
        rf_im = rf_im(1 + floor(-sz_diff/2):end - ceil(-sz_diff/2), :, :, :);
    end
        
        
    input_im{1} = I;
    input_im{2} = -input_im{1};
    
    time_of_simulation = ceil((size(I,2) + size(rf_im,2)) / abs(v) + size(neuron_irf,1));
    
    conv_sig = zeros(time_of_simulation, size(rf_im, 3), 2);
            
    
    
    for l = 1:2
        
        for n = 1:size(rf_im, 3)                        
            dx = -Inf;
            for t = 1:time_of_simulation
                
                old_dx = dx;
                if v > 0
                    dx = ceil(t*v);
                else
                    dx = floor(t*v) + size(rf_im, 2) + size(I,2);
                end

                if old_dx ~= dx
                    
                    I_start = max(dx - size(rf_im, 2) + 1, 1);
                    I_end = min(dx, size(I,2));

                    rf_end = min(size(rf_im, 2) + size(I,2) - dx, size(rf_im,2));
                    rf_start = max(1, size(rf_im, 2) - dx + 1);

                    conv_sig(t, n, l) = sum(sum( ...
                        rf_im(:, rf_start:rf_end, n, l) .* ...
                        input_im{l}(:, I_start:I_end) ));                    
                else
                    conv_sig(t, n, l) = conv_sig(t-1, n, l);
                end

                %
            end
            
            temp_sig = conv(conv_sig(:, n, l), neuron_irf(:, n));
            conv_sig(:, n, l) = temp_sig(1:size(conv_sig,1));            
        end
    end
    
    conv_sig(conv_sig < 0) = 0;
    for n = 1:size(neuron_irf,3)
        conv_sig(:,n,l) = conv_sig(:,n,l) * in_sign(n);
    end
    
    out_sig = sum(conv_sig(:,:), 2);
end
