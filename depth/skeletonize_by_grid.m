function [nodes edges] = skeletonize_by_grid(im, chunk_size)
    
    MAX_CUBE_PTS = 32;
    MAX_POINTS = 2048;

    imSz = size(im);
    
    num_chunks = ceil(imSz./chunk_size);
    
    is_explored = false(num_chunks);

    start_point = find(im(:), 1, 'first'); 
    [start_point(1) start_point(2) start_point(3)] = ind2sub(imSz, start_point);
    
    
    nodes = zeros(MAX_POINTS,3);
    edges = zeros(MAX_POINTS,2);
    
    to_explore = ceil(start_point./chunk_size);
    f_nums = [2 3; 1 3; 1 2; 2 3; 1 3; 1 2];
    f_d = [1; 2; 3; 1; 2; 3];
    f_nhood = [-eye(3); eye(3)];
    
    total_node_k = 0;
    total_edge_k = 0;
    
    while ~isempty(to_explore);
        tic
        disp(['exploring cube ' num2str(to_explore(1,:))]);
        
        cube_begin = (to_explore(1,:)-1).*chunk_size + 1;
        cube_end = (to_explore(1,:)).*chunk_size + 1; %1 more for overlap
        cube_end = min([cube_end; imSz]);
        
        cubeSz = cube_end-cube_begin+1;
        
        cube = im(cube_begin(1):cube_end(1), ...
            cube_begin(2):cube_end(2), ...
            cube_begin(3):cube_end(3));
            
        
        cube_conns = bwconncomp(cube,26);
        
        cube_pts = zeros(MAX_CUBE_PTS,3);
        cube_edges = zeros(MAX_CUBE_PTS,2);
        cube_k = 0;
        edge_k = 0;
        
        f_bounds = [1 1 1 cubeSz];
        
        
        
        for ob = 1:cube_conns.NumObjects
            
            [x(:,1) x(:,2) x(:,3)] = ind2sub(cubeSz,cube_conns.PixelIdxList{ob});          
            
            x_is_bound = [x == 1, x(:,1)==f_bounds(4), x(:,2)==f_bounds(5), x(:,3)==f_bounds(6)];    
            
            ob_k = 0;
            
            
            %commented below: treat each face independently
%             for f = 1:6;
%                 if any(x_is_bound(:,f))
%                     f_im = zeros(cubeSz(f_nums(f,:)));
%                     f_inds = sub2ind(cubeSz(f_nums(f,:)), x(x_is_bound(:,f),f_nums(f,1)), x(x_is_bound(:,f),f_nums(f,2)));
% 
%                     f_im(f_inds) = true;
%                     f_conns = bwconncomp(f_im,8);
% 
%                     for f_ob = 1:f_conns.NumObjects
%                         [y z] = ind2sub(cubeSz(f_nums(f,:)), f_conns.PixelIdxList{f_ob});
%                         ob_k = ob_k+1;
%                         cube_pts(cube_k+ob_k,f_nums(f,1)) = mean(y);
%                         cube_pts(cube_k+ob_k,f_nums(f,2)) = mean(z);
%                         cube_pts(cube_k+ob_k,f_d(f)) = f_bounds(f);
%                     end
%                     
%                     to_explore = [to_explore; to_explore(1,:) + f_nhood(f,:)];
%                 end                                                                
%             end
            
            %new plan is to do cc on the entire shell
            
            bx = x(any(x_is_bound,2),:);
            f_im = false(cubeSz);
            f_im(sub2ind(cubeSz,bx(:,1),bx(:,2),bx(:,3))) = true;
            f_conns = bwconncomp(f_im,26);
            for f_ob = 1:f_conns.NumObjects
                [p q r] = ind2sub(cubeSz, f_conns.PixelIdxList{f_ob});
                ob_k = ob_k+1;
                cube_pts(cube_k+ob_k,1) = mean(p);
                cube_pts(cube_k+ob_k,2) = mean(q);
                cube_pts(cube_k+ob_k,3) = mean(r);
            end
            
            for f = 1:6
                if any(x_is_bound(:,f))
                    to_explore = [to_explore; to_explore(1,:) + f_nhood(f,:)];
                end
            end
            


            if ob_k == 1 %if only 1 point in object, also get furthest
                dist_from_point = (x(:,1)-cube_pts(1,1)).^2 + ...
                    (x(:,2)-cube_pts(1,2)).^2 + ...
                    (x(:,3)-cube_pts(1,3)).^2;
                [~, max_ind] = max(dist_from_point);
                ob_k = 2;
                cube_pts(cube_k+ob_k,:) = x(max_ind,:);
                
            end
            
            if ob_k == 2 %make edge
                edge_k = edge_k+1;
                cube_edges(edge_k,:) = [cube_k+ob_k-1 cube_k+ob_k];
            elseif ob_k >= 3 %make new centroid point, make edge from it to others
                cube_pts(cube_k+ob_k+1,:) = ...
                    mean(cube_pts(cube_k+(1:ob_k),:));
                
                cube_edges(edge_k+(1:ob_k),1) = cube_k+(1:ob_k);
                cube_edges(edge_k+(1:ob_k),2) = cube_k+ob_k+1;
                
                ob_k = ob_k+1;
                edge_k = edge_k+ob_k-1;                                
            end
            
            cube_k = cube_k+ob_k;
            
            x = [];
        end
        
        %prepare local lists for global lists
        cube_pts = cube_pts(1:cube_k,:); 
        cube_edges = cube_edges(1:edge_k,:);
        
        for d = 1:3
            cube_pts(:,d) = cube_pts(:,d) + cube_begin(d) - 1;
        end
        cube_edges = cube_edges + total_node_k;
        
        edges(total_edge_k+(1:edge_k),:) = cube_edges;
        nodes(total_node_k+(1:cube_k),:) = cube_pts;
        
        total_edge_k = total_edge_k+edge_k;
        total_node_k = total_node_k+cube_k;
        
        
        %handle exploration
        is_explored(to_explore(1,1), to_explore(1,2), to_explore(1,3)) = true;
        for d = size(to_explore,1):-1:1
            if any(to_explore(d,:)==0) || any(to_explore(d,:)>num_chunks) || ...
                    is_explored(to_explore(d,1), to_explore(d,2), to_explore(d,3))
                to_explore(d,:) = [];
            end
        end
        disp(['nodes found: ' num2str(cube_k) '; edges found: ' num2str(edge_k) '; time elapsed: ' num2str(toc)]);
        
    end
    
    nodes = nodes(1:total_node_k,:);
    edges = edges(1:total_edge_k,:);
    
    
end
        

    