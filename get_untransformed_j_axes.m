function get_untransformed_j_axes
    C = get_constants;
    load('~/stratification/T.mat');
    for n = C.type.j
        
        load(['~/stratification/surface_points/cell_' num2str(n) '_surface.mat']);
        cell_dat = cell_data(n);
        soma = cell_dat.get_midpoint(true);
    
        
        




%         surface_points(C.f(surface_points(1,:)) > 80,:) = [];
        
        surface_points = double(surface_points);
        
        for d = 1:3
            surface_points(:,d) = surface_points(:,d)*C.res(d);
        end
        
        surface_points = surface_points*T.global_Q';
        
        min_p = min(surface_points(:,1));
        
        soma = mean(surface_points(surface_points(:,1)<min_p+400,:));
        
        surface_points = surface_points(:,2:3);
        
        
        
        for d = 1:2
            surface_points(:,d) = (surface_points(:,d)-soma(d+1));
        end

        distal_p = get_distal_loc(n);
            
            daxis = distal_p(2:3) - soma(2:3);
            daxis = daxis/norm(daxis);
            theta = atan2(daxis(2),daxis(1));
            step = pi/8;
            step_dir = 1;

            num_p = size(surface_points,1);

            while 1
                orth_d = [-daxis(2) daxis(1)];
                num_pos = sum(surface_points(:,1)*orth_d(1) + surface_points(:,2)*orth_d(2)>0);
%                     disp(num_pos/num_p)
                if abs(num_pos/num_p - .5) < .001
                    break
                elseif num_pos/num_p < .5
                    %too few in pos direction, sweep ccw
                    if step_dir == 1
                        step = step/2;
                        step_dir = -1;
                    end
                elseif num_pos/num_p > .5
                    if step_dir == -1
                        step = step/2;
                        step_dir = 1;
                    end
                end

                theta = theta + step_dir*step;
                daxis = [cos(theta) sin(theta)];

            end
            
        
            
            
        dist_point = soma(2:3) + 100000*daxis;
        
        dist_point = [soma(1) dist_point];
        dist_point = dist_point*T.global_Q;
        soma = soma * T.global_Q;
        
        soma = soma./C.res;
        dist_point = dist_point ./ C.res;
        
%         dist_point = dist_point/C.res(2);
        
        
%         plot([soma(2); dist_point(1)], [soma(3); dist_point(2)],'Color', [1 0 0],'LineWidth', 3);
%         soma = 1;
        soma = round(soma); 
        dist_point = round(dist_point);
        disp([num2str(n) '    (' num2str(soma(1)) ', ' num2str(soma(2)) ', ' num2str(soma(3)) ')   (' num2str(dist_point(1)) ', ' num2str(dist_point(2)) ', ' num2str(dist_point(3)) ')'])
        
        
        
    end