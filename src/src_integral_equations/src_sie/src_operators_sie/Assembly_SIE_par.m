function [Z,Z_copper_loss] = Assembly_SIE_par(SCOIL,emc,inputs) 

    % -------------------------------------------------------------------------
    %            Info 
    % -------------------------------------------------------------------------

    index = SCOIL.index;
    etod = SCOIL.etod;
    node = SCOIL.node;
    elem = SCOIL.elem;
    Index_elem = SCOIL.index_elem;

    Nvars  = max(index);

    % -------------------------------------------------------------------------
    %     Order of Gauss Quadrature Integration 
    % -------------------------------------------------------------------------

    N_NS         = inputs.N_NS; 
    N_ST_psi     = inputs.N_ST_psi;
    N_EA_theta   = inputs.N_EA_theta; 
    N_EA_psi     = inputs.N_EA_psi;
    N_VA_theta_p = inputs.N_VA_theta_p;
    N_VA_theta_q = inputs.N_VA_theta_q;
    N_VA_psi     = inputs.N_VA_psi;

    GL_order.ST = N_ST_psi;
    GL_order.EA = [N_EA_theta N_EA_psi];
    GL_order.VA = [N_VA_theta_p N_VA_theta_q N_VA_psi];
    GL_order.NS  = N_NS;

    % -------------------------------------------------------------------------
    %            Start the assembly procedure 
    % -------------------------------------------------------------------------

    Z = zeros(Nvars,Nvars);

    % -------------------------------------------------------------------------
    %             NS interactions
    % -------------------------------------------------------------------------

    [Z] = assembly_ns_par(index,etod,node,elem,Z,GL_order,Index_elem,emc.k0);

    % -------------------------------------------------------------------------
    %             EA interactions 
    % -------------------------------------------------------------------------

    [Z] = assembly_ea_par(index,etod,node,elem,Z,GL_order,Index_elem,emc.k0);

    % -------------------------------------------------------------------------
    %             VA interactions 
    % -------------------------------------------------------------------------

    [Z] = assembly_va_par(index,etod,node,elem,Z,GL_order,Index_elem,emc.k0);

    Z = Z + Z.'; 

    % -------------------------------------------------------------------------
    %             ST interactions 
    % -------------------------------------------------------------------------

    [Z,Z_copper_loss] = assembly_st_par(index,etod,node,elem,Z,GL_order,Index_elem,emc);

    Z = (emc.eta0 / (4 * pi)) * Z;
    Z_copper_loss = (emc.eta0 / (4 * pi)) * Z_copper_loss; 

end

