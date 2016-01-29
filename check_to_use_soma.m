function use_soma = check_to_use_soma(cell_num)
    C = get_constants;
    
    bip = [C.type.BC1 C.type.BC2 C.type.BC3a C.type.BC3b C.type.BC4 C.type.on_bc];
    
    if any(bip==cell_num)
        use_soma = false;
    else
        use_soma = true;
    end
        
end