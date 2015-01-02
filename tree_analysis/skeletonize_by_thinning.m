function skeletonize_by_thinning(cell_num, varargin)

    p = inputParser;    
    p.addRequired('cell_num', @isnumeric);
    p.addOptional('chunk_size', [300 300 300], @(x) isnumeric(x) && length(x) == 3);
    p.addOptional('expected_max_radius', 5, @isnumeric);    
    
    p.parse(cell_num, varargin{:});    
    s = p.Results;

    
    C = get_constants;
    
    cell_dat = cell_data(cell_num);        
    bb = cell_dat.get_volume_bounds;
    
    num_chunks = ceil((bb(:,2)' - bb(:,1)' - s.expected_max_radius)./(s.chunk_size-2*s.expected_max_radius));
    
    tot_chunks = prod(num_chunks);
    
    [chunk_loc(:,1), chunk_loc(:,2), chunk_loc(:,3)] = ind2sub(num_chunks, 1:tot_chunks);             
    
    nhood1 = [eye(3); -eye(3)];
    
    [nhood2(:,1), nhood2(:,2), nhood2(:,3)] = ind2sub([3 3 3], 1:prod([3 3 3]));
    nhood2 = nhood2-2;
    nhood2 = nhood2(any(nhood2,2),:);
    
    is_simple = false(2^size(nhood2,1),1);
    
%     tic
%     for n = 1:size(is_simple,1)-1
%         
%         
%         K = false(3,3,3);
%         for k = 1:size(nhood2,1)
%             K(2+nhood2(k,1), 2+nhood2(k,2), 2+nhood2(k,3)) = mod(floor((n-1)/(2^(k-1))),2)==1;
%         end
%         
%         segs = bwconncomp(K,size(nhood2,1));
%         is_simple(n) = segs.NumObjects == 1;
%         
%     end
%     toc
%     save('./skele_nhood.mat','is_simple');
    
    load('./skele_nhood.mat');
    
    
    for n = 1:tot_chunks
        
        cstart = bb(:,1)' + (chunk_loc(n,:)-1).*(s.chunk_size-2*s.expected_max_radius);
        cend = cstart + s.chunk_size - 1;
        
        p = cell_dat.get_volume([cstart; cend]');
        
        
        if ~isempty(p)

            for d = 1:3
                p(:,d) = p(:,d) - cstart(d) + 1;
            end

            data_cube = true(s.chunk_size+2);
            data_cube(2:end-1, 2:end-1, 2:end-1) = false;

            
            for k = 1:size(p,1)
                data_cube(p(k,1)+1, p(k,2)+1, p(k,3)+1) = true;
            end

            disp(['working on chunk ' num2str(n) ' of ' num2str(tot_chunks)]);
            tic
            for t = 1:s.expected_max_radius-1
                
                is_surface = false(s.chunk_size);

                for k = 1:size(nhood1,1);
                    is_surface( ...
                        data_cube((2:end-1) + nhood1(k,1), ...
                        (2:end-1) + nhood1(k,2), ...
                        (2:end-1) + nhood1(k,3)) ~= ...
                        data_cube(2:end-1, 2:end-1, 2:end-1)) = true;
                end
                is_surface(~data_cube(2:end-1, 2:end-1, 2:end-1)) = false;
                
                surface_ind = find(is_surface(:));
                new_p = zeros(length(surface_ind),3);
                [new_p(:,1) new_p(:,2) new_p(:,3)] = ind2sub(size(is_surface), surface_ind);
                
                new_p = new_p+1; %put into data_cube coords
                
                for l = 1:size(new_p,1)
                    nhood_type = 0;
                    for k = 1:size(nhood2,1);
                        nhood_type = nhood_type + 2^(k-1) * ...
                            data_cube(new_p(l,1) + nhood2(k,1), ...
                            new_p(l,2) + nhood2(k,2), ...
                            new_p(l,3) + nhood2(k,3));
                    end
                    if nhood_type~=0 && is_simple(nhood_type)
                        data_cube(new_p(l,1), new_p(l,2), new_p(l,3)) = false;
                    end
                    
                end
                
                
                nhood_type = ones(s.chunk_size);
                nhood_type(nhood_type==0) = size(is_simple,1);
                new_cube = data_cube(2:end-1, 2:end-1, 2:end-1) & ~is_simple(nhood_type);
                
                

                sum(data_cube(:))
                data_cube(2:end-1, 2:end-1, 2:end-1) = new_cube;
                sum(new_cube(:))
                
                

            end
            toc
            
            cube_p = find(data_cube(:));
            [new_p(:,1) new_p(:,2) new_p(:,3)] = ind2sub(size(data_cube), cube_p);

            for d = 1:3
                new_p(:,d) = new_p(:,d)+bb(d,1)-1+t;
            end
            skele_p = [skele_p; new_p];

        
        end
        
    end
            
    