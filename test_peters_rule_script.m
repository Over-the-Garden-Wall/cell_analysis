function test_peters_rule_script

    C = get_constants;
    figure; hold all
    
    bip_types = {'t1' 't2' 't3a' 't3b' 't4'};
    
    load(C.conn_loc);
    
    reliable_threshold = 2;
    
    bin_size = 5*10^7;
    
    for k = 1:5
        
        ol = find_symmetric_overlap(C.type.off_sac, C.type.(bip_types{k}));
        
        is_reliable = false(size(ol));
        obs_conts = zeros(size(ol));
        
        for n = 1:length(C.type.off_sac)
            R = pick_conns(conns, C.type.off_sac(n), C.type.(bip_types{k}));
    
            for m = 1:length(R)
                if size(R{m},2) >= reliable_threshold
                    obs_conts(n,m) = sum(double(R{m}(3,:)));
                    is_reliable(n,m) = true;
                end
            end
        end
        
        ol = (ol(is_reliable));
        obs_conts = (obs_conts(is_reliable));
                
        
        x = 0:bin_size:ceil(max(ol));
        y = zeros(size(x));
        for n = 1:length(x);
            v = ceil(ol/bin_size)*bin_size==x(n);
            if any(v)
                y(n) = mean(obs_conts(v));
            end
        end
        
        plot(x,y);
        
%         scatter(ol, obs_conts, '*');
        
        
    end
end