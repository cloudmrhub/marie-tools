function[MREDM] = ie_solver_svie_mrgf(MREDM,U_hat_inv,S_hat_inv,V_hat_inv,X)

    dims = MREDM.dimensions;
    res  = dims.res;
    ql   = dims.ql;
    l    = dims.l;
    iG   = reshape(repelem([1;12]/res^3,[1,l-1],3),1,[]);

    MREDM.solver.rhs = rhs_assembly(MREDM,0);

    % Wire-Coil-Shield
    if ~isempty(MREDM.WIE.coil) && ~isempty(MREDM.SIE.coil) && ~isempty(MREDM.SIE.shield)
        Zss  = MREDM.SIE.Z_shield;
        Zcc  = MREDM.SIE.Z;
        Zww  = MREDM.WIE.Z;
        Zsw  = MREDM.operators.Zsw.U * MREDM.operators.Zsw.V';
        Zsc  = MREDM.operators.Zsc.U * MREDM.operators.Zsc.V';
        Zwc  = MREDM.operators.Zcw.U * MREDM.operators.Zcw.V';
        Zall = [Zss   Zsw   Zsc;
                Zsw.' Zww   Zwc;
                Zsc.' Zwc.' Zcc];

        ZbcN = [MREDM.operators.Zbs_N;
                MREDM.operators.Zbw_N;
                MREDM.operators.Zbc_N];
    % Wire-Coil   
    elseif ~isempty(MREDM.WIE.coil) && ~isempty(MREDM.SIE.coil) && isempty(MREDM.SIE.shield)
        Zcc  = MREDM.SIE.Z;
        Zww  = MREDM.WIE.Z;
        Zwc  = MREDM.operators.Zcw.U * MREDM.operators.Zcw.V';
        Zall = [Zww   Zwc;
                Zwc.' Zcc];

        ZbcN = [MREDM.operators.Zbw_N;
                MREDM.operators.Zbc_N];
    % Coil-Shield
    elseif isempty(MREDM.WIE.coil) && ~isempty(MREDM.SIE.coil) && ~isempty(MREDM.SIE.shield)
        Zss  = MREDM.SIE.Z_shield;
        Zcc  = MREDM.SIE.Z;
        Zsc  = MREDM.operators.Zsc.U * MREDM.operators.Zsc.V';
        Zall = [Zss   Zsc;
                Zsc.' Zcc];

        ZbcN = [MREDM.operators.Zbs_N;
                MREDM.operators.Zbc_N];
    % Coil
    elseif isempty(MREDM.WIE.coil) && ~isempty(MREDM.SIE.coil) && isempty(MREDM.SIE.shield)
        Zcc  = MREDM.SIE.Z;
        Zall = Zcc;

        ZbcN = MREDM.operators.Zbc_N;
    % Wire-Shield
    elseif ~isempty(MREDM.WIE.coil) && isempty(MREDM.SIE.coil) && ~isempty(MREDM.SIE.shield)
        Zss  = MREDM.SIE.Z_shield;
        Zww  = MREDM.WIE.Z;
        Zsw  = MREDM.operators.Zsw.U * MREDM.operators.Zsw.V';
        Zall = [Zss   Zsw;
                Zsw.' Zww];

        ZbcN = [MREDM.operators.Zbs_N;
                MREDM.operators.Zbw_N];
    % Wire
    elseif ~isempty(MREDM.WIE.coil) && isempty(MREDM.SIE.coil) && isempty(MREDM.SIE.shield)
        Zww  = MREDM.WIE.Z;
        Zall = Zww;

        ZbcN = MREDM.operators.Zbw_N;
    end

    % MRGF Solution
    tic
    Zss = Zall + (ZbcN).' * reshape(iG.*reshape(U_hat_inv*(S_hat_inv*(V_hat_inv'*ZbcN)),[],ql,N_coil),[],N_coil);
    Jcb = Zss\MREDM.solver.rhs;
    a   = X*(reshape(iG.*reshape(ZbcN,[],ql,N_coil),[],N_coil)*Jcb);
    fprintf('\tMRGF Solution: Time: %s \n',hh_mm_ss(toc));
    
    % Store
    MREDM.fields.Jcb = Jcb;
    MREDM.fields.a   = a;

    fprintf('\tComputing network parameters ...\n');
    MREDM = np_compute(MREDM);

end
