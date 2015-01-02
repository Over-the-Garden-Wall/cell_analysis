function [x, y, does_intersect] = line_intersect(line_xy_a, line_xy_b, intsct_xy_a, intsct_xy_b)

    if line_xy_a(1) > line_xy_b(1)
        temp = line_xy_b;
        line_xy_b = line_xy_a;
        line_xy_a = temp;
    end
    
    needs_swapping = intsct_xy_a(:,1) > intsct_xy_b(:,1);
    temp = intsct_xy_b(needs_swapping,:);
    intsct_xy_b(needs_swapping,:) = intsct_xy_a(needs_swapping,:);
    intsct_xy_a(needs_swapping,:) = temp;
    

    line_slope = (line_xy_b(2)-line_xy_a(2))/(line_xy_b(1)-line_xy_a(1));
    int_slopes = (intsct_xy_b(:,2)-intsct_xy_a(:,2))./(intsct_xy_b(:,1)-intsct_xy_a(:,1));
    
    int_x = ((intsct_xy_a(:,2) - intsct_xy_a(:,1).*int_slopes) - ...
        (line_xy_a(2) - line_xy_a(1)*line_slope))./(line_slope-int_slopes);
    
    int_y = line_slope*(int_x-line_xy_a(1)) + line_xy_a(2);
    
    does_intersect = int_x >= intsct_xy_a(:,1) & int_x <= intsct_xy_b(:,1) & ...
        int_x >= line_xy_a(1) & int_x <= line_xy_b(1); 
    
    x = int_x(does_intersect);
    y = int_y(does_intersect);
    
end
    