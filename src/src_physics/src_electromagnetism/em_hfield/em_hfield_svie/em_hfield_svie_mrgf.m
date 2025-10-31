function[h_field] = em_hfield_svie_mrgf(MREDM,Ub)

    mu0 = MREDM.emc.mu0;
    
    % Compute Total Field
    dims         = MREDM.dimensions;
    n1           = dims.n1;
    n2           = dims.n2;
    n3           = dims.n3;
    ql           = dims.ql;
    l            = dims.l;
    idx_solution = dims.idx_solution;

    for i = 1:2
        if ~isempty(MREDM.fields.JcbTx) && i == 1
            ac = MREDM.fields.aTx;
        elseif ~isempty(MREDM.fields.JcbRx) && i == 2        
            ac = MREDM.fields.aRx;
        else
            continue
        end
        A_ports = size(ac,2);
    
        h_b                 = zeros(n1*n2*n3*ql,A_ports);
        h_b(idx_solution,:) = Ub*ac;
        h_b                 = reshape(h_b,n1,n2,n3,ql,A_ports);
        b1m                 = mu0 * (squeeze(h_b(:,:,:,1,:)) - 1i*squeeze(h_b(:,:,:,1+l,:)));
        b1p                 = mu0 * (squeeze(h_b(:,:,:,1,:)) + 1i*squeeze(h_b(:,:,:,1+l,:)));
  
        % Save to struct and return
        if ~isempty(MREDM.fields.JcbTx) && i == 1
            h_field.Tx.h_b = h_b;
            h_field.Tx.b1m = b1m;
            h_field.Tx.b1p = b1p;
        elseif ~isempty(MREDM.fields.JcbRx) && i == 2      
            h_field.Rx.h_b = h_b;
            h_field.Rx.b1m = b1m;
            h_field.Rx.b1p = b1p;
        end
    end

end
