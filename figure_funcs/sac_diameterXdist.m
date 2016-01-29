%70183
%70023
C = get_constants;


cns = C.type.on_sac;
% cns(1:find(cns==70183,1,'first')) = [];

% cns = C.type.off_sac;
% cns(1:find(cns==70023,1,'first')) = [];


x = 0:5:200;

y = zeros(size(x));
yn = zeros(size(x));
figure; hold all

for cn = 1:5; %length(cns)
    
    y = zeros(size(x));
    yn = zeros(size(x));

    c = cns(cn);
    
    s = load([C.skele_dir 's' num2str(c) '.mat']);
    
    c_d = cell_data(c);
    soma_loc = c_d.get_midpoint(true);
    
    z = C.f(s.nodes(:,1));
    s.nodes(z < 20 | z > 75,:) = [];
    s.node_diameter(z < 20 | z > 75) = [];
    
    d = sqrt( (s.nodes(:,2) - soma_loc(2)).^2 + (s.nodes(:,3) - soma_loc(3)).^2 );
    
    d = floor(d/1000 / (x(2)-x(1)) )+1;
    
    for n = 1:length(d)
        y(d(n)) = y(d(n)) + s.node_diameter(n);
        yn(d(n)) = yn(d(n)) + 1;
    end
    plot(x, y./yn);
end

% plot(x, y./yn);
        
    
    
    