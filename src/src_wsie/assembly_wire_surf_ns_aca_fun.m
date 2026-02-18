function [U,V] = assembly_wire_surf_ns_aca_fun(first_p, second_p, third_p, index_surface, etod_surface, node_surface, elem_surface, Np_1D, Np_2D, k0, tol)
    
    % Constants and quadratures
    first_node             = [3 1 2];
    second_node            = [2 3 1];
    [Np,wt,Z_sie,z1,z2,z3] = gauss_2d_sie(Np_2D);
    [wr,Z_w]               = gauss_1d(Np_1D);
    W_t                    = wr*wt;
    Zi_j                   = zeros(Np_1D,Np,3);
    A_o                    = 1i*k0*ones(Np_1D,Np);
    F_o                    = 1/(1i*k0)*4*ones(Np_1D,Np);
    Ne_surface             = size(elem_surface, 2);
    Ne_wire                = length(third_p);
    
    % Precompute the Zi_j integrals (for all combinations of quadrature points)
    % Dimensions: Np1D x Np2D x 3.
    for jj=1:3
        Zi_j(:,:,jj) = (Z_w + 1)*Z_sie(:,jj).';
    end

    % Entries for all wire points and Precompute for columns
    ro_x_left  = zeros(Np_1D, Ne_wire);
    ro_y_left  = zeros(Np_1D, Ne_wire);
    ro_z_left  = zeros(Np_1D, Ne_wire);
    ro_x_right = zeros(Np_1D, Ne_wire);
    ro_y_right = zeros(Np_1D, Ne_wire);
    ro_z_right = zeros(Np_1D, Ne_wire);
    fn_left    = zeros(1, Ne_wire);
    fn_right   = zeros(1, Ne_wire);
    tlp_left   = zeros(3, Ne_wire);
    tlp_right  = zeros(3, Ne_wire);
    Dl         = vecnorm(second_p(1:Ne_wire,:)-first_p(1:Ne_wire,:),2,2);
    Dr         = vecnorm(third_p(1:Ne_wire,:)-second_p(1:Ne_wire,:),2,2);
    xi_shifted = (Z_w + 1); 
    for n = 1:Ne_wire
        % Nodes
        p1 = first_p(n,:);
        p2 = second_p(n,:);
        p3 = third_p(n,:);
        % Position in Segment
        lp_left  = (Dl(n)/2) * xi_shifted;     
        lp_right = (Dr(n)/2) * xi_shifted;   
        % Real Position
        ro_x_left(:,n)  = p1(1) + (lp_left / Dl(n)) * (p2(1) - p1(1));
        ro_y_left(:,n)  = p1(2) + (lp_left / Dl(n)) * (p2(2) - p1(2));
        ro_z_left(:,n)  = p1(3) + (lp_left / Dl(n)) * (p2(3) - p1(3));
        ro_x_right(:,n) = p2(1) + (lp_right / Dr(n)) * (p3(1) - p2(1));
        ro_y_right(:,n) = p2(2) + (lp_right / Dr(n)) * (p3(2) - p2(2));
        ro_z_right(:,n) = p2(3) + (lp_right / Dr(n)) * (p3(3) - p2(3));
        % Define constant tangent vectors (the segments are straight)
        tlp_left(:,n)  = [(p2(1)-p1(1)), (p2(2)-p1(2)), (p2(3)-p1(3))] / Dl(n);
        tlp_right(:,n) = [(p3(1)-p2(1)), (p3(2)-p2(2)), (p3(3)-p2(3))] / Dr(n);
        % Testing Functions         
        fn_left(1,n)  = (Dl(n)/2) / Dl(n);            
        fn_right(1,n) = (Dr(n) - (Dr(n)/2)) / Dr(n);  
    end
    roX_left  = reshape(ro_x_left,  [Np_1D, 1, Ne_wire, 1]); 
    roY_left  = reshape(ro_y_left,  [Np_1D, 1, Ne_wire, 1]);
    roZ_left  = reshape(ro_z_left,  [Np_1D, 1, Ne_wire, 1]);
    roX_right = reshape(ro_x_right, [Np_1D, 1, Ne_wire, 1]); 
    roY_right = reshape(ro_y_right, [Np_1D, 1, Ne_wire, 1]);
    roZ_right = reshape(ro_z_right, [Np_1D, 1, Ne_wire, 1]);

    % Entries for coil points
    atsd         = abs(etod_surface);
    node_basis_1 = elem_surface(1,:);
    node_basis_2 = elem_surface(2,:);
    node_basis_3 = elem_surface(3,:);
    rs_1         = node_surface(:,node_basis_1);
    rs_2         = node_surface(:,node_basis_2);
    rs_3         = node_surface(:,node_basis_3);
    L_src        = zeros(3,3,Ne_surface);
    L_src(:,1,:) = rs_2-rs_3;
    L_src(:,2,:) = rs_3-rs_1;
    L_src(:,3,:) = rs_1-rs_2;
    rs_x         = z1.'*rs_1(1,:)+z2.'*rs_2(1,:)+z3.'*rs_3(1,:);
    rs_y         = z1.'*rs_1(2,:)+z2.'*rs_2(2,:)+z3.'*rs_3(2,:);
    rs_z         = z1.'*rs_1(3,:)+z2.'*rs_2(3,:)+z3.'*rs_3(3,:);
    rsX          = reshape(rs_x, [1, Np, 1, Ne_surface]); 
    rsY          = reshape(rs_y, [1, Np, 1, Ne_surface]);
    rsZ          = reshape(rs_z, [1, Np, 1, Ne_surface]);

    % Precompute for Rows
    idx_surface = [];
    j_surface   = [];
    ssj         = [];
    j_1         = [];
    j_2         = [];
    je_idx      = [];
    for je = 1:Ne_surface
        as = squeeze(atsd(:,je));
        for j=1:3
            if index_surface(as(j)) ~= 0
                j_surface   = [j_surface;j];
                idx_surface = [idx_surface;index_surface(as(j))];
                ssj      = [ssj;sign(etod_surface(j, je))];
                j_1      = [j_1;second_node(j)];
                j_2      = [j_2;first_node(j)];
                je_idx   = [je_idx; je];
            end
        end
    end
    global_idx1   = (je_idx - 1) * 3 + j_surface; 
    global_idx_j1 = (je_idx - 1) * 3 + j_1;
    global_idx_j2 = (je_idx - 1) * 3 + j_2;
    L_src_flat    = reshape(L_src, 3, []); 
    L_src_vecs    = L_src_flat(:, global_idx1);  
    L_j           = sqrt(sum(L_src_vecs.^2, 1)).'; 
    L_src_vecs_j1 = L_src_flat(:, global_idx_j1); 
    L_src_vecs_j2 = L_src_flat(:, global_idx_j2); 

    % Initialize ACA
    I = 1;
	U = [];
	V = [];

    ZR_row = zeros(1,max(index_surface));
    ZR_col = zeros(Ne_wire,1);
	
	%Ri = zeros(M, 1);
	SF2 = 0;

    for k=1:min(Ne_wire,Ne_surface)
        % Compute rows
        Ri = assembly_wire_surf_ns_row(Dl, Dr, fn_left, fn_right, tlp_left, tlp_right, ro_x_left, ro_y_left, ro_z_left, ro_x_right, ro_y_right, ro_z_right, rsX, rsY, rsZ, Np_1D, Np, W_t, Zi_j, A_o, F_o, k0, ZR_row, idx_surface, je_idx, L_j, L_src_vecs_j1, L_src_vecs_j2, j_1, j_2, ssj, I);
        if k > 1
            Ri = Ri - U(I,:)*V';
        end
        [~,J] = max(abs(Ri));
        Vk = (Ri / Ri(J))';
        
        % Compute columns
        Uk = assembly_wire_surf_ns_column(Ne_surface, index_surface, etod_surface, atsd, L_src, roX_left, roY_left, roZ_left, roX_right, roY_right, roZ_right, rs_x, rs_y, rs_z, Np_1D, Np, W_t, reshape(Zi_j,Np_1D,Np,1,3), A_o, F_o, first_node, second_node, k0, ZR_col, Ne_wire, Dl, Dr, fn_left, fn_right, tlp_left, tlp_right, J);
        
        if k > 1
            Uk = Uk - U*(V(J,:))';
        end

        % Compute Approximation
        nVk = norm(Vk);
        nUk = norm(Uk);
        SF2 = SF2 + (nVk*nUk)^2; 
        if k > 1
            SF2 = SF2 + 2*sum(real((U'*Uk).*(V'*Vk))); 
        end

        % Evaluate if we are done
        U = [U Uk];
        V = [V Vk];
        if nUk*nVk <= tol*sqrt(SF2)
            break;
        end
        Uk = abs(Uk);
        Uk(I) = 0;
        [~,I] = max(Uk);
        
    end
    
end
