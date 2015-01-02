function sanity_check

    d = [.1 .2 .3 .4 .5];
    
%     d = [.4 .5];
    
    exp = 2*(d'*d);
    obs = zeros(5);
    
    for m = 1:5
        
        for n = 1:5
            
            
            for k = 1:10
                r = rand(100,100,100);
                im = -(r <= d(m)) + (r >= 1-d(n));                                
                
                obs(m,n) = obs(m,n) + ...
                    mean(mean(mean(-1 == im(1:end-1,:,:).*im(2:end,:,:)))) + ...
                    mean(mean(mean(-1 == im(:,1:end-1,:).*im(:,2:end,:)))) + ...
                    mean(mean(mean(-1 == im(:,:,1:end-1).*im(:,:,2:end))));
            end
            obs(m,n) = obs(m,n)/3/10;
        end
    end
    
    figure;scatter(exp(:), obs(:));
    
    figure;scatter(sqrt(exp(:)), obs(:));
end