function ps = apply_transform(T,ps)
    BATCH_SIZE = 1000000;
    
    for n = 1:size(ps,2)
        ps(:,n) = ps(:,n)*T.res(n);
    end

    ps = ps*T.global_Q';
        
    if strcmp(T.transform, 'original')
        for n = 1:size(ps,1)
            p = ps(n,:);

            dist = sqrt((p(2)-T.nodes(:,1)).^2 + ...
                (p(3)-T.nodes(:,2)).^2);

            valid_inds = dist<= T.reach;
        %     valid_inds = (dist==min(dist));

            if ~any(valid_inds)
                [~, d_ind] = min(dist);
                q = T.quaternion(d_ind,:);
                translatn = T.translation(d_ind,:);
            else
                weight = T.reach - dist(valid_inds);
                weight = weight/sum(weight);                

                q = sum(T.quaternion(valid_inds,:).*(weight*ones(1,4)),1);
                translatn = sum(T.translation(valid_inds,:).*(weight*ones(1,3)), 1);
            end

        %     Q = angles2Q(theta, phi, psi);
            Q = quaternion2Q(q);

            p = p*Q' + translatn;
            ps(n,:) = p;
        end
            
    elseif strcmp(T.transform, 'orthogonal')

        for n = 1:size(ps,1)
            p = ps(n,:);
            dist = sqrt((p(2)-T.nodes(:,1)).^2 + ...
            (p(3)-T.nodes(:,2)).^2);

            valid_inds = dist<= T.reach;
        %     valid_inds = (dist==min(dist));

            if ~any(valid_inds)
                [~, d_ind] = min(dist);
                translatn = T.translation(d_ind,:);
            else
                weight = T.reach - dist(valid_inds);
                weight = weight/sum(weight);                

                translatn = sum(T.translation(valid_inds,:).*(weight*ones(1,3)), 1);
            end

            p = p + translatn;
            ps(n,:) = p;          
        end
    elseif strcmp(T.transform, 'orthostretch');
        num_total_points = size(ps,1);
        ps2 = ps;
        for n = 1:ceil(num_total_points/BATCH_SIZE);
            
            tps = ps(((n-1)*BATCH_SIZE + 1):min(n*BATCH_SIZE, num_total_points),:);
            num_points = size(tps,1);
            
            num_nodes = size(T.nodes,1);
            dist = sqrt((tps(:,2)*ones(1,num_nodes) - ones(num_points,1)*T.nodes(:,1)').^2 + ...
                (tps(:,3)*ones(1,num_nodes) - ones(num_points,1)*T.nodes(:,2)').^2);
            valid_inds = dist<= T.reach;
            has_valid = any(valid_inds,2);

            w = T.reach-dist;
            w(~has_valid,:) = w(~has_valid,:)+T.reach;
            w(w<0) = 0;


            t = (ones(num_points, 1)*T.off_chat') .* w;
            t = sum(t,2) ./ sum(w,2);
            s = (ones(num_points, 1)*(T.off_chat-T.on_chat)') .* w;
            s = sum(s,2) ./ sum(w,2);

        
            ps(((n-1)*BATCH_SIZE + 1):min(n*BATCH_SIZE, num_total_points),1) = ...
                (tps(:,1) - t)./s * T.mean_band_distance;
        end
        
%         for n = 1:size(ps,1)
%             p = ps(n,:);
%             dist = sqrt((p(2)-T.nodes(:,1)).^2 + ...
%             (p(3)-T.nodes(:,2)).^2);
% 
%             valid_inds = dist<= T.reach;
%         %     valid_inds = (dist==min(dist));
%             if ~any(valid_inds)
%                 [~, d_ind] = min(dist);                
%                 translatn = T.off_chat(d_ind,:);
%                 stretch = T.on_chat(d_ind,:)-T.off_chat(d_ind,:);                
%             else                
%                 weight = T.reach - dist(valid_inds);
%                 weight = weight/sum(weight);                
% 
%                 translatn = sum(T.off_chat(valid_inds).*weight);
%                 stretch = sum((T.off_chat(valid_inds)-T.on_chat(valid_inds)).*weight);   
%             end
% 
%             p(1) = (p(1) - translatn)/stretch * T.mean_band_distance;
%             ps(n,:) = p;          
%         end
    
    end
end