function annos = parse_omni_annotation(fn)

    store_len = 10;
    

%example below
% - id: 82
%   enabled: true
%   value: {coord: [2826, 8510, 6268], comment: 70119, color: [255, 0, 0], size: 3}

    fid = fopen(fn);
    store = [];
    annos = [];
    field_name = [];
    while 1
        new_char = fread(fid,1,'*char');
        if isempty(new_char)
            break
        end
        
        if length(store) < store_len
            store = [store new_char];
        else
            store = [store(2:end) new_char];
        end
        
        
        if new_char == ':'
            name_begin = find(store==' ' | store=='{',1,'last')+1;
            if isempty(name_begin)
                field_name = store(1:end-1);
            else
                field_name = store(name_begin:end-1);
            end
            new_char = fread(fid,1,'*char');
            collecting_array = false;
            store = [];
        elseif new_char == '['
            collecting_array = true;
            store = [];
            vals = [];
        elseif ~isempty(field_name) && (new_char == ',' || new_char == '}' || new_char == 10)
            if ~collecting_array
                annos = append_value(annos, field_name, str2double(store(1:end-1)));
            else
                vals(end+1) = str2double(store(1:end-1));
                new_char = fread(fid,1,'*char');                
                store = [];
            end
        elseif new_char == ']'
            vals(end+1) = str2double(store(1:end-1));
            annos = append_value(annos, field_name, vals);
            
        end
    end
        
    
    
    
    fclose(fid);
end

function a = append_value(a,fn, v)
    if isempty(fn)
        a=a;
    elseif isfield(a,fn)
        a.(fn)(end+1,:) = v;
    else
        a.(fn) = v;
    end
end
