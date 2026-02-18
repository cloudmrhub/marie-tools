function [M_new,F_new,N_new] = map_to_local(M,F,N)

    all_ids = [F(:);M(:);N(:)];
    
    [sorted_ids, ~, sorted_idx] = unique(all_ids);
    [~, rank_map] = sort(sorted_ids);  % indices into sorted order
    ranking = zeros(1, numel(sorted_ids));
    ranking(rank_map) = 1:numel(sorted_ids);
    
    % Apply mapping
    mapped = ranking(sorted_idx);
    F_new = mapped(1:length(F));
    M_new = mapped(length(F)+1:length(M)+length(F));
    N_new = mapped(length(M)+length(F)+1:end);

end