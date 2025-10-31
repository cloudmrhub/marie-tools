function [portentry, loop_start, loop_end, F_point, S_point, T_point] = Mesh_Wire(fname)
    
    fid = fopen(fname,'r');
    C = textscan(fid,'%s','Delimiter','\n');
    fclose(fid);
    lines = C{1};
    
    % Extract $Nodes block
    iN1 = find(strcmp(lines,'$Nodes'))+1;
    nNodes = str2double(lines{iN1});
    nodeLines = lines(iN1+1:iN1+nNodes);
    nodes = zeros(nNodes,3);
    for i=1:nNodes
        parts = sscanf(nodeLines{i},'%d %f %f %f');
        nodes(i,:) = parts(2:4).';
    end
    
    % Extract $Elements block
    iE1 = find(strcmp(lines,'$Elements'))+1;
    nElems = str2double(lines{iE1});
    elemLines = lines(iE1+1:iE1+nElems);
    
    ports = [];
    F_point = []; S_point = [];
    L1 = []; L2 = [];
    
    for i=1:nElems
        parts = sscanf(elemLines{i},'%d %d %d %d %d %d %d');
        etype = parts(2);
        if etype==15  % point element
            nid = parts(end);
            ports(end+1) = nid; %#ok<AGROW>
        elseif etype==1  % line element
            n1 = parts(end-1);
            n2 = parts(end);
            L1(end+1) = n1; %#ok<AGROW>
            L2(end+1) = n2; %#ok<AGROW>
            F_point(end+1,:) = nodes(n1,:); %#ok<AGROW>
            S_point(end+1,:) = nodes(n2,:); %#ok<AGROW>
        end
    end
    
    % mimic old variables
    portentry  = ports(:);
    if ~isempty(S_point)
        T_point = S_point(2:end,:);
    else
        T_point = [];
    end

    [loop_start, loop_end] = find_starts_ends(L1,L2);

end

function [loop_start, loop_end] = find_starts_ends(L1, L2)
    % L1, L2: [N x 1] segment connectivity (node IDs)
    % Returns segment indices: loop_start, loop_end
    
    N = numel(L1);
    starts = [];
    ends   = [];
    
    allNodes = [L1(:); L2(:)];
    uniqNodes = unique(allNodes);
    
    % ---------- Loops ----------
    for n = uniqNodes.'
        idxL1 = find(L1 == n);
        idxL2 = find(L2 == n);
        if ~isempty(idxL1) && ~isempty(idxL2)
            i = idxL1(1);                                 % earliest start at node n
            pCand = idxL2(abs(idxL2 - i) > 1);            % enforce separation > 1
            if ~isempty(pCand)
                % prefer the first L2 occurrence after the start; fallback to first with separation
                p = pCand(find(pCand > i, 1, 'first')); 
                if isempty(p), p = pCand(1); end
                starts(end+1) = i; %#ok<AGROW>
                ends(end+1)   = p; %#ok<AGROW>
            end
        end
    end
    
    % ---------- Dipoles (open chains) ----------
    % nodes that appear exactly once across L1+L2 (loose ends)
    maxNode = max(allNodes);
    counts  = accumarray(allNodes, 1, [maxNode, 1]);
    looseNodes = find(counts == 1);
    
    if ~isempty(looseNodes)
        % build interleaved sequence: L1(1),L2(1), L1(2),L2(2), ...
        seqNodes  = zeros(2*N,1);
        seqRow    = zeros(2*N,1);  % segment index for that entry
        seqIsL1   = false(2*N,1);  % true if from L1, false if from L2
        for i = 1:N
            k = 2*i - 1;
            seqNodes(k) = L1(i);   seqRow(k) = i; seqIsL1(k) = true;
            seqNodes(k+1) = L2(i); seqRow(k+1) = i; % L2
        end
    
        % positions of each loose node in interleaved sequence
        pos = zeros(numel(looseNodes),1);
        for t = 1:numel(looseNodes)
            pos(t) = find(seqNodes == looseNodes(t), 1, 'first');
        end
        [~, ordLoose] = sort(pos);
        looseOrdered = looseNodes(ordLoose);
    
        % pair consecutive loose nodes as (start,end) in order
        for k = 1:2:numel(looseOrdered)
            if k+1 > numel(looseOrdered), break; end
            nStart = looseOrdered(k);
            nEnd   = looseOrdered(k+1);
    
            % start index: prefer where it appears in L1; fallback to its L2 row
            iStart = find(L1 == nStart, 1, 'first');
            if isempty(iStart), iStart = find(L2 == nStart, 1, 'first'); end
    
            % end index: prefer where it appears in L2; fallback to its L1 row
            iEnd = find(L2 == nEnd, 1, 'last');
            if isempty(iEnd), iEnd = find(L1 == nEnd, 1, 'last'); end
    
            % append
            starts(end+1) = iStart; %#ok<AGROW>
            ends(end+1)   = iEnd;   %#ok<AGROW>
        end
    end
    
    % ---------- Clean up: dedupe & order by appearance ----------
    % If a node produced multiple (start,end) pairs, keep earliest start
    [starts, ia] = unique(starts, 'stable');
    ends = ends(ia);
    
    [~, order] = sort(starts, 'ascend');
    loop_start = starts(order)';
    loop_end   = ends(order)';
end
















