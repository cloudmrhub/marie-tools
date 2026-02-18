function [M_new,N_new] = map_to_local_2(M,N)

    all_ids = [M(:);N(:)];
    
    [sorted_ids, ~, sorted_idx] = unique(all_ids);
    [~, rank_map] = sort(sorted_ids);  % indices into sorted order
    ranking = zeros(1, numel(sorted_ids));
    ranking(rank_map) = 1:numel(sorted_ids);
    
    % Apply mapping
    mapped = ranking(sorted_idx);
    M_new = mapped(1:length(M));
    N_new = mapped(length(M)+1:end);

end