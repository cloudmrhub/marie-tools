function [ZR] = assembly_ns_par(index,etod,node,elem,ZR,GL_order,Index_elem,ko)
%%

%
A_o    = 1i*ko;
F_o    = 1/(1i*ko);
%
first_node  = [3 1 2];
second_node = [2 3 1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                Order of Gauss Quadrature Integration                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Np_2D  = GL_order.NS;   
%%%%%%%%%%%%%%% Gauss quadrature rule for non-singular triangles %%%%%%%%%%
[Np,wt,Z,~,~,~] = gauss_2d_sie(Np_2D);

% generate a weighted combination of quadratures
Z_Wt = diag(wt') * Z;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% sort list by observer triangles
Ne = size(elem,2);
NS_indices = Index_elem.NS(:, 1); 
NS_values = Index_elem.NS(:, 2); 
NS_cell = accumarray(NS_indices, NS_values, [Ne, 1], @(x) {x});
standard_empty = zeros(0, 1);
for i = 1:numel(NS_cell)
    if isempty(NS_cell{i})
        NS_cell{i} = standard_empty;
    end
end

% n_NS_elem = size(Index_elem.NS,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             Assembly                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% time_ns = tic;
%%%%%%%%%%%%%%%%%%%%%%%%% Main body of assembly %%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%  loop through all triangles
%%%%%%%%%%%%%%%%%%%%%%%%% Non-Singular Terms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get total number of unique triangles
N_tr  = size(NS_cell,1);

% allocate memory for auxiliary matrices
r      = zeros(N_tr * Np, 3);
l      = zeros(3 * N_tr, 3);
Ln_ij  = zeros(N_tr, 3);

% generate 3x3 grid with indices to unfold inner-most loops
[i_inner, j_inner] = meshgrid(1:3, 1:3);
i_inner = reshape(i_inner,9,1);
j_inner = reshape(j_inner,9,1);

% index offset 
index_offset = reshape(3 * ones(9,1) * (0:N_tr-1), 9*N_tr,1);

% generate indices for unfolded loop over src triangles
i_1 = repmat(second_node(i_inner)', N_tr, 1);
i_2 = repmat(first_node(i_inner)', N_tr, 1);
j_1 = repmat(second_node(j_inner)', N_tr, 1); 
j_2 = repmat(first_node(j_inner)', N_tr, 1);

% indices for L_ij_so; Note, that j's and i's have to be flipped, since src
% index corresponds to the row and an observer index to the column
j1_i1 = sub2ind([3*N_tr,3], j_1 + index_offset, i_1);
j2_i1 = sub2ind([3*N_tr,3], j_2 + index_offset, i_1);
j1_i2 = sub2ind([3*N_tr,3], j_1 + index_offset, i_2);
j2_i2 = sub2ind([3*N_tr,3], j_2 + index_offset, i_2);

% "global" indices for L and W structures
L_glob_idx = [j1_i1 j2_i1 j1_i2 j2_i2];
W_glob_idx = [j2_i2 j1_i2 j2_i1 j1_i1];

% generate indices to fetch edge orientations (signs)
obs_sign_idx = repmat(i_inner, N_tr,1);

% index shifts due to reduction in number of interacting src triangles
L_shift   = 3 * [i_1 i_1 i_2 i_2];
W_shift   = 3 * [i_2 i_2 i_1 i_1];

% vectors for reshape and inner products
LW_mult   = [1 -1 -1 1]';
nine_ones = ones(9,1);
col_plus  = ones(N_tr * Np,1);
row_plus  = ones(1,Np); 

% precompute auxiliary matrices
for i = 1:N_tr
    n1 = elem(1,i);
    n2 = elem(2,i);
    n3 = elem(3,i);
    %
    r_1 = node(:,n1);
    r_2 = node(:,n2);
    r_3 = node(:,n3);
    
    ind_range_Np = Np*(i-1)+1:Np*i;
    ind_range_3  = 3*(i-1)+1:3*i;
    
    r(ind_range_Np,:) = Z * [r_1, r_2, r_3]';
    l(ind_range_3,:) = [r_2-r_3, r_3-r_1, r_1-r_2]';
    Ln_ij(i,:) = sqrt(diag(l(ind_range_3,:) * l(ind_range_3,:)'));
end

x = r(:,1); y = r(:,2); z = r(:,3);
x_obs = x.'; y_obs = y.'; z_obs = z.'; 

a = [x.^2, x, col_plus,...
     y.^2, y, col_plus,...
     z.^2, z, col_plus];
     
% loop across observer triangles
for io = 1:N_tr

    % get src triangles for a given observer triangle
    src_triangles = nonzeros(NS_cell{io}); 
    
    Ntr_i = N_tr - io;
    lsi   = Ntr_i * Np;
    N_src = length(src_triangles);

    src_range = io*Np+1: N_tr*Np;

    % physical edges for current observer triangle
    index_ao = index(repmat(abs(etod(i_inner,io)),N_src,1));
    index_as = index(reshape(abs(etod(j_inner,src_triangles)),9*N_src,1));
    
    % define range of src triangles within a vector of all triangles
    obs_ind_range = (io-1) * Np + 1: io * Np;

    b = [row_plus; -2 * x_obs(obs_ind_range); x_obs(obs_ind_range).^2; ...
         row_plus; -2 * y_obs(obs_ind_range); y_obs(obs_ind_range).^2; ...
         row_plus; -2 * z_obs(obs_ind_range); z_obs(obs_ind_range).^2];
     
    c = a(src_range,:);
    
    % compute Dinstances and weighted Green's Functions 
    G = sqrt(c * b);

    G     = reshape(exp(-1i*ko*G)./G,lsi,Np); % Calculate Ns*Np by Np Greens func.
    wGw   = 4 * F_o * reshape(nine_ones * (wt * reshape(G*wt',Np,Ntr_i)),[],1);  % A row vector of Ns scalars.
    zwGzw = A_o * reshape(Z_Wt'*reshape(G*Z_Wt, Np, 3*Ntr_i), 9*Ntr_i,1);
    
    L_ij  = reshape(l(3*io+1:end,:) * l(3*(io-1)+1:3*io,:)', 9*Ntr_i,1);

    L_tr_idx = L_glob_idx - io * L_shift;
    W_tr_idx = W_glob_idx - io * W_shift;

    L_loc_idx = reshape(L_tr_idx(9*io+1:end,:), 36 * Ntr_i,1);
    W_loc_idx = reshape(W_tr_idx(9*io+1:end,:), 36 * Ntr_i,1);
    

    L_ij12   = reshape(L_ij(L_loc_idx), 9*Ntr_i,4);
    W_GR_A_o = reshape(zwGzw(W_loc_idx), 9*Ntr_i,4);   
    
    % get vectors of 9*Ntr_i edge lenghts
    L_i = repmat(Ln_ij(io, i_inner)',Ntr_i,1);
    L_j = reshape(Ln_ij(io+1:end, j_inner)', 9*Ntr_i,1);
    
    % get orientations of edges in observer and src triangles
    soi = sign(etod(obs_sign_idx(1:end-9 * io), io));
    ssj = reshape(sign(etod(j_inner,io+1:end)), 9*Ntr_i,1);
    
    Z_io = reshape(soi .* ssj .* L_i .* L_j .*((L_ij12 .* W_GR_A_o) * LW_mult + wGw), 9, Ntr_i);
    
    % get indices of src triangles with respect to io offset
    src_idx = src_triangles - io;

    inner_mask = find(index_ao .* index_as);
    
    asidx = index_as(inner_mask);
    aoidx = index_ao(inner_mask);
    
    idx_1d = sub2ind(size(ZR), asidx, aoidx);
    [ind_1d,~,unique_mask] = unique(idx_1d);
    
    Z_val = reshape(Z_io(:,src_idx), [],1);
    z_vec =  Z_val(inner_mask);
    v = accumarray(unique_mask, z_vec);
    
    ZR(ind_1d) = ZR(ind_1d) + v;

end

clear Z_io Z_val z_vec asidx aoidx;