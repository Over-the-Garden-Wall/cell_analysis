function [classification, all_types] = classify_cells_by_full_strat(unknown_cells, run_mode)
    
if ~exist('unknown_cells','var')
    unknown_cells = [];
end

if ~exist('run_mode', 'var') || isempty(run_mode)
    if isempty(unknown_cells)
        run_mode = 'validate';
    else
        run_mode = 'classify';
    end
end

% num_cells = length(cell_nums);
C = get_constants;

types.unknown = unknown_cells;
num_unknown = length(unknown_cells);

%types as defined by jinseop
types.t1 = [60008, 60019, 60026, 60027, 60032, 60052, 60055, 60078, 60079, 60099, 60105, 60109, 60110, 60111, 60114, 60118, 60129, 60132, 60139, 60142, 60147, 60150, 60158, 60161, 60162, 60164, 60170, 60177, 60184, 60187, 60188, 60189, 60194, 60195, 60196, 60203, 60204, 60213, 60216, 60218, 60227];
types.t2 = [60001, 60002, 60003, 60004, 60006, 60009, 60010, 60011, 60012, 60013, 60021, 60022, 60023, 60025, 60029, 60037, 60038, 60039, 60040, 60041, 60042, 60043, 60046, 60080, 60097, 60101, 60102, 60103, 60104, 60112, 60120, 60124, 60130, 60133, 60135, 60138, 60140, 60141, 60149, 60157, 60159, 60167, 60182, 60208, 60209, 60210, 60219, 60221, 60223, 60224, 60226, 60198];
types.t3a = [60028, 60030, 60048, 60049, 60050, 60056, 60059, 60066, 60068, 60072, 60075, 60076, 60085, 60088, 60092, 60093, 60108, 60115, 60123, 60134, 60136, 60137, 60145, 60146, 60166, 60172, 60176, 60181];
types.t3b = [60024, 60045, 60053, 60057, 60063, 60069, 60073, 60074, 60077, 60082, 60084, 60087, 60090, 60094, 60096, 60098, 60106, 60107, 60122, 60131, 60148, 60151, 60153, 60154, 60156, 60165, 60173, 60178, 60179, 60190, 60201, 60211, 60215, 60228, 60081];
types.t4 = [60047, 60058, 60067, 60083, 60086, 60089, 60091, 60095, 60100, 60116, 60117, 60119, 60121, 60125, 60127, 60143, 60144, 60160, 60163, 60168, 60171, 60174, 60175, 60183, 60185, 60186, 60191, 60197, 60199, 60202, 60206, 60220, 60222];
types.t5w = [60054, 60071, 60155, 60366, 60371, 60395, 60399, 60403, 60405, 60408, 60411, 60436, 60440, 60445, 60452, 60465, 60481, 60502, 60503, 60507, 60527, 60532, 60533, 60535, 60540];
types.t5l = [60035, 60020, 60033, 60360, 60374, 60380, 60383, 60386, 60388, 60389, 60404, 60410, 60414, 60415, 60421, 60439, 60442, 60450, 60458, 60460, 60462, 60478, 60488, 60490, 60497, 60498, 60504, 60505, 60510, 60513, 60514, 60519, 60522, 60523, 60528, 60541, 60542];
types.t5h = [60018, 60061, 60064, 60387, 60390, 60406, 60420, 60432, 60434, 60447, 60448, 60453, 60459, 60466, 60469, 60473, 60484, 60491, 60495, 60534, 60538, 60543];
types.t6 = [60036, 60060, 60356, 60358, 60363, 60370, 60375, 60381, 60382, 60394, 60398, 60401, 60407, 60409, 60416, 60423, 60428, 60431, 60441, 60443, 60444, 60451, 60463, 60464, 60475, 60476, 60477, 60483, 60486, 60487, 60489, 60496, 60499, 60512, 60516, 60526, 60530, 60536, 60537, 60548, 60549, 60554, 60614];
types.t7 = [60016, 60051, 60354, 60354, 60361, 60373, 60376, 60377, 60393, 60396, 60397, 60412, 60418, 60429, 60435, 60437, 60446, 60454, 60470, 60480, 60485, 60492, 60508, 60509, 60525, 60529, 60531, 60546, 60553, 60562];
types.t89 = [60368, 60402, 60417, 60422, 60433, 60438, 60457, 60461, 60482, 60500, 60561, 60578];
types.tR = [60017, 60031, 60357, 60359, 60365, 60369, 60372, 60378, 60384, 60392, 60400, 60456, 60468, 60471, 60472, 60474, 60479, 60506, 60518, 60520, 60521, 60544, 60551, 60555, 60557, 60558, 60560, 60563, 60564, 60565, 60566, 60567, 60568, 60569, 60570, 60571, 60572, 60573, 60574, 60575, 60576, 60577, 60579, 60580, 60581, 60582, 60583, 60584, 60585, 60586, 60587, 60588, 60589, 60590, 60591, 60592, 60593, 60594, 60595, 60596, 60597, 60598, 60599, 60600, 60601, 60602, 60603, 60604, 60605, 60606, 60607, 60609, 60610, 60613];
types.tX = [60355, 60379, 60413, 60430, 60449, 60455, 60493, 60501, 60517, 60539, 60547, 60550];


type_names = fieldnames(types);
num_types = length(type_names);

all_gt_nums = [];
all_types = [];
for t = 1:num_types
    all_gt_nums = [all_gt_nums types.(type_names{t})];
    all_types(end+1:length(all_gt_nums)) = 2^(t-1);    
end
num_total_cells = length(all_gt_nums);


data_matrix = zeros(103, num_total_cells);

%construct stratification profiles
for t = 1:num_total_cells
    c = all_gt_nums(t);
    c_d = cell_data(c);
    
    p = c_d.get_surface;
    d = round(C.f(p(:,1)));
    
    is_valid = d >= 1 & d <= 100;
    
    p = p(is_valid,:);
    d = d(is_valid);
    
    coverage_p = round(p(:,2:3));
    data_matrix(101, t) = size(unique(coverage_p, 'rows'), 1);
    data_matrix(102, t) = c_d.V;
    
    h = cell(1,2);
    [h{:}] = poly2cw(c_d.hull_2d(:,1), c_d.hull_2d(:,2));
    data_matrix(103, t) = polyarea(h{1}, h{2});
    
    s = c_d.stratification(1-C.strat_x(1):end);
    if length(s) > 100
        s = s(1:100);
    end
    data_matrix(1:length(s),t) = s/sum(s);
        
end

for n = 1:3
    data_matrix(100+n,:) = data_matrix(100+n,:)/max(data_matrix(100+n,:));
end

data_matrix = data_matrix';

classification = zeros(size(all_types));
    
if strcmp(run_mode, 'validate')
    
    for k = 1:num_total_cells        
        former_type = all_types(k);
        all_types(k) = 1;
        
        classifier = generate_classifier(data_matrix(all_types~=1, :), all_types(all_types~=1));
        classification(k) = run_classifier(classifier, data_matrix(k,:));
        all_types(k) = former_type;
    end

    
else
    
    classifier = generate_classifier(data_matrix(all_types~=1, :), all_types(all_types~=1));
    classification = run_classifier(classifier, data_matrix);
            
end

end

%pairwise svm
function svms = generate_classifier(M, lbl)

    lbl_nums = unique(lbl);
    lbl_nums(lbl_nums==1) = [];
    
    svms = cell(length(lbl_nums));
    
    for k = 1:length(lbl_nums)
        for j = k+1:length(lbl_nums)
            is_valid = lbl == lbl_nums(k) | lbl == lbl_nums(j);
            svms{k, j} = svmtrain(M(is_valid,:), 2*(lbl(is_valid) == lbl_nums(k)) - 1, 'autoscale', 'false');
        end
    end

end

function clss = run_classifier(svms, M)

    svmout = zeros([size(M,1), size(svms)]);
    
    for k = 1:size(svms,1)
        for j = k+1:size(svms,2)
            svmout(:,k, j) = svmclassify(svms{k,j}, M);
            %             svmout(:,j, k) = (svms{k,j}.SupportVectors * M')'*svms{k,j}.Alpha(:) + svms{k,j}.Bias;
        end
    end
    
    
    votes = zeros(size(M,1), size(svms,1));
    
    for k = 1:size(svms,1)
        votes(:,k) = sum(svmout(:,k, :), 3) - sum(svmout(:,:, k), 2);
    end
    
    clss = zeros(size(svmout,1),1);
    for k = 1:size(svmout,1)
        [num_votes, max_ind] = max(votes(k,:)); 
        if sum(votes(k,:)==num_votes) > 1
            candidates = votes(k,:) == num_votes;
            sub_svmout = squeeze(svmout(k, candidates, candidates));
            subvotes = zeros(size(sub_svmout,1),1);
            for n = 1:size(sub_svmout,1)
                subvotes(n) = sum(sub_svmout(k, :), 2) - sum(sub_svmout(:, k), 1);
            end
            [d, max_ind] = max(subvotes(k,:)); 
            tmp = find(candidates, max_ind, 'first');
            clss(k) = tmp(end);
        else
            clss(k) = max_ind;
        end
    end
        
    
end



% %winner takes all svm
% function model = generate_classifier(M, lbl)
% 
%     lbl_nums = unique(lbl);
%     lbl_nums(lbl_nums==1) = [];
%     is_valid = lbl~=1;
%     
%     for k = 1:length(lbl_nums)
%         model{k} = svmtrain(M(is_valid,:), lbl(is_valid) == lbl_nums(k), 'autoscale', 'false');
%     end
% 
% end
% 
% function clss = run_classifier(svms, M)
% 
%     clss = zeros(size(M,1), length(svms));
%     
%     for k = 1:length(svms)
%         clss(:,k) = M * svms{k}.SupportVectors' * svms{k}.Alpha(:) + svms{k}.Bias;
%     end
%     
%     [clss_val, clss] = min(clss, [], 2); 
%     
% end

    