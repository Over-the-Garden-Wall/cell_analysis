function [subskeles node_correspondance] = split_skele(skele, split_nodes)

    G = skeleton2graph(skele);
    G(split_nodes,:) = 0;
    G(:,split_nodes) = 0;

    [S, C] = graphconncomp(G);

    subskeles = cell(S,1);
    node_correspondance = cell(S,1);
    for n = 1:S
        node_subset = find(C==n);
        node_correspondance{n} = node_subset;
        sub_G = G(node_subset,node_subset);

        subskeles{n}.edges = graph2edgeList(sub_G);
        subskeles{n}.nodes = skele.nodes(node_subset,:);

    end

    for n = S:-1:1
        if isempty(subskeles{n}.edges)
            subskeles(n) = [];
            node_correspondance(n) = [];
        end
    end

end
