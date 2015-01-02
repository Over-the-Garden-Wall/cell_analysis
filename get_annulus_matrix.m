function A = get_annulus_matrix(p1, p2, num_annuls1, num_annuls2)
    %A is a matrix of the areas of annuli overlap, where Aij represents the
    %total area of the overlap of annulus i around p1 with annulus j around
    %p2
    
    if length(p1) == 3
        p1 = p1(2:3);
    end
    if length(p2) == 3
        p2 = p2(2:3);
    end
    
    
    A = zeros(num_annuls1, num_annuls2);
    
    d = sqrt(sum((p2 - p1).^2))/1000;
    
    for m = 1:num_annuls1
        for n = 1:num_annuls2
            
                A(m,n) = get_circle_intersect(d, m, n) - ...
                    get_circle_intersect(d, m-1, n) - ...
                    get_circle_intersect(d, m, n-1) + ...
                    get_circle_intersect(d, m-1, n-1);
            
            
        end
    end
    
    
end

function A = get_circle_intersect(d, r1, r2) 
    if d > r1 + r2
        A = 0;
    elseif d <= r2 - r1
        A = r1^2*pi;
    elseif d <= r1 - r2
        A = r2^2*pi;
    else
        CBA = acos((r2^2 + d^2 - r1^2)/(2*r2*d));
        CBD = 2*CBA;
        
        CAB = acos((r1^2 + d^2 - r2^2)/(2*r1*d));
        CAD = 2*CAB;
        
         A = r2^2*(CBD - sin(CBD))/2 + ...
             r1^2*(CAD - sin(CAD))/2;
         
%          if imag(A)~= 0
%              A = A;
%          end
         
    end
        
end
    



    