function [list_a, list_b] = flip_list_items(list_a, list_b, items_to_flip)
    if size(list_a,2) > 1
        list_a = list_a';
    end
    if size(list_b,2) > 1
        list_b = list_b';
    end
    if size(items_to_flip,1) > 1
        items_to_flip = items_to_flip';
    end
    
    
    for t = items_to_flip
        if any(t == list_a)
            list_b = [list_b; t];
            list_a(t==list_a) = [];
        else
            list_a = [list_a; t];
            list_b(t==list_b) = [];
        end
    end
end
            