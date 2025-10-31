function [Z] = Assembly_SIE_block_par(SCOIL,SHIELD,emc,inputs) 

    index_coil        = SCOIL.index;
    etod_coil         = SCOIL.etod;
    node_coil         = SCOIL.node;
    elem_coil         = SCOIL.elem;
    index_shield      = SHIELD.index;
    etod_shield       = SHIELD.etod;
    node_shield       = SHIELD.node;
    elem_shield       = SHIELD.elem;
    
    [U,Z.V] = assembly_surf_ns_aca_fun(index_shield,etod_shield,node_shield,elem_shield,index_coil,etod_coil,node_coil,elem_coil,inputs.N_NS,emc.k0,inputs.tol_ACA);        
    Z.U = -(emc.eta0 / (4 * pi)) * U;

end