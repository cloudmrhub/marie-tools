function [permutedElem] = Mesh_Permute(elem,e)
    permutedElem = elem;
    % Loop through each line in e
    for i = 1:size(e,2)
        node1 = e(1, i);
        node2 = e(2, i);
        % Find the two triangles containing both node1 and node2
        triangles = [];
        for j = 1:size(elem,2)
            nodes = elem(1:3, j);
            if ismember(node1, nodes) && ismember(node2, nodes)
                triangles = [triangles, j]; %#ok<AGROW>
            end
        end
        
        if length(triangles) ~= 2
            error('Expected to find exactly two triangles for each line.');
        end
        
        % Get the nodes of the two triangles
        tri1 = elem(1:3, triangles(1));
        
        % Permute the triangles to align the connection
        % Find the order of nodes in tri1 that matches the connection node1 -> node2
        idx1 = find(tri1 == node1);
        idx2 = find(tri1 == node2);
        if mod(idx2 - idx1, 3) ~= 1
            % Swap the triangles in permutedElem if the order does not align
            temp = permutedElem(:, triangles(1));
            permutedElem(:, triangles(1)) = permutedElem(:, triangles(2));
            permutedElem(:, triangles(2)) = temp;
        end
    end
end