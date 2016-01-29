C = get_constants;
% 
% type_is_predefined = true;
% output_is_coverage = true;
is_on_dsgc = false;
% use_face_only = true;

skip_factor = 10;

if ~type_is_predefined
    if is_on_dsgc
        dsgcs = C.type.on_dsgc;
    else
        dsgcs = C.type.oodsgc;
    end
end

num_gc = length(dsgcs);

% if type_is_predefined
%     pref_axes = ones(num_gc,2);
% elseif is_on_dsgc

    tic
    pref_axes = zeros(num_gc,2);
    for n = 1:num_gc
        pref_axes(n,:) = dsgc_preferred_direction(dsgcs(n))';
    end
    toc
    pref_axes
% else
% 
%     pref_axes = [-0.9962   -0.0872
%    -0.4067    0.9135
%    -0.9135    0.4067
%    -0.0349   -0.9994
%     0.9744    0.2250
%     0.9976    0.0698
%    -0.9816    0.1908
%     0.3090   -0.9511
%     0.9998    0.0175
%     0.9994    0.0349
%    -0.8290   -0.5592
%    -0.9744   -0.2250
%     0.9877    0.1564];
% end 

if use_face_only
    load(C.raw_conn_loc);
end

dsgc_group = zeros(num_gc,1);
for n = 1:num_gc
    [dummy, coord] = max(abs(pref_axes(n,:)));
    dsgc_group(n) = (round(pref_axes(n,coord))+1) / 2 + coord*2 - 1;
end

group_dir = zeros(4,1);
for n = 1:4
    ax = mean(pref_axes(dsgc_group==n,:),1);
    
    group_dir(n) = atan2(ax(2), ax(1));
end

BC_types = {'t1', 't2', 't3a', 't3b', 't4', 't5w', 't5l', 't5h', 't6', 't7'};
type_layer = [1 1 1 1 1 2 2 2 2 2];

BCs = []; BC_typenum = []; BC_layer = [];
for type_k = 1:length(BC_types);
    BCs = [BCs C.type.(BC_types{type_k})];
    BC_typenum(end+1:length(BCs)) = type_k;
    BC_layer(end+1:length(BCs)) = type_layer(type_k);
end
num_bc = length(BCs);

bc_dist = zeros(num_gc, num_bc);
bc_contact = zeros(num_gc, num_bc);
bc_overlap = zeros(num_gc, num_bc);

l_p = cell(2,1);
for dk = 1:num_gc;
    d = dsgcs(dk);
        
    
    c_d = cell_data(d);
    
    if use_face_only
        is_me = conns(1,:) == d;
        d_contacts = double(conns(2:5,is_me));
        
        is_me = conns(2,:) == d;
        d_contacts = [d_contacts, double(conns([1 3:5],is_me))];
        
        d_contacts(2,:) = sum(d_contacts(2:4,:));
    else        
        d_contacts = double(c_d.contacts);
    end
    
    p = c_d.get_surface;
    pdepth = C.f(p(:,1));
    
    l_p{1} = p(pdepth < 46.5,2:3);
    l_p{2} = p(pdepth >= 46.5 & pdepth < 70, 2:3);
    
    l_mid{1} = mean(l_p{1});
    l_mid{2} = mean(l_p{2});
    
    total_p = size(l_p{1},1) + size(l_p{2},1);
    l_porportion = [size(l_p{1},1) / total_p, size(l_p{2},1) / total_p];
    
    l_p{1} = l_p{1}(ceil(1:skip_factor:end),:);
    l_p{2} = l_p{2}(ceil(1:skip_factor:end),:);
    
    gc_theta = group_dir(dsgc_group(dk));
    gc_axis = [cos(gc_theta), sin(gc_theta)];
    
    tic
    for bk = 1:num_bc;
        b = BCs(bk);
        bc_data = cell_data(b);
        
        l = BC_layer(bk);
        
        bc_hull = [];
        [bc_hull(:,1), bc_hull(:,2)] = poly2cw(bc_data.hull_2d(:,1), bc_data.hull_2d(:,2));
        
        
        bc_contact(dk, bk) = sum(d_contacts(2, d_contacts(1,:)==b));
        
        bc_overlap(dk, bk) = sum(inpolygon(l_p{l}(:,1), l_p{l}(:,2), ...
            bc_hull(:,1), bc_hull(:,2))) * skip_factor;
        
        bc_mid = bc_data.get_midpoint;
        
        bc_dist(dk, bk) = sum((bc_mid(2:3) - l_mid{l}) .* gc_axis);
        
        if output_is_coverage
            %correction so that contact/overlap is representative of total
            %dsgc coverage
            bc_overlap(dk, bk) = bc_overlap(dk, bk) / l_porportion(l);
            
        end
        
    end
    toc
    
end


%connection density per cell
type_contact = zeros(num_gc, length(BC_types));
type_overlap = zeros(num_gc, length(BC_types));
for n = 1:length(BC_types)
    type_contact(:,n) = sum(bc_contact(:, BC_typenum==n), 2);
    type_overlap(:,n) = sum(bc_overlap(:, BC_typenum==n), 2);
end
figure; plot(type_contact./type_overlap);

contxtypexloc = zeros(length(BC_types), 4, 2);
overlapxtypexloc = zeros(length(BC_types), 4, 2);

for b = 1:length(BC_types)
    btype = BC_typenum(b);
    for gct = 1:4
        my_gcs = dsgc_group == gct;
        
        sub_bc_contact = bc_contact(my_gcs, BC_typenum==b);
        sub_bc_overlap = bc_overlap(my_gcs, BC_typenum==b);
        
        sub_dist = bc_dist(my_gcs, BC_typenum==b);
        
        contxtypexloc(b, gct, 1) = sum(sum(sub_bc_contact(sub_dist>0)));
        contxtypexloc(b, gct, 2) = sum(sum(sub_bc_contact(sub_dist<0)));
        
        overlapxtypexloc(b, gct, 1) = sum(sum(sub_bc_overlap(sub_dist>0)));
        overlapxtypexloc(b, gct, 2) = sum(sum(sub_bc_overlap(sub_dist<0)));
        
        
    end
end
        
figure; plot(sum(contxtypexloc,3) ./ sum(overlapxtypexloc,3));    
       
figure; plot((contxtypexloc(:,:,1) - contxtypexloc(:,:,2)) ./ ...
    (overlapxtypexloc(:,:,1) + overlapxtypexloc(:,:,2)));