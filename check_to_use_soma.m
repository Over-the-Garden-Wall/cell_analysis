function use_soma = check_to_use_soma(cell_num)
    C = get_constants;
    
    bip = [C.type.t1 C.type.t2 C.type.t3a C.type.t3b C.type.t4];
    
    if any(bip==cell_num)
        use_soma = false;
    else
        use_soma = true;
    end
        
end