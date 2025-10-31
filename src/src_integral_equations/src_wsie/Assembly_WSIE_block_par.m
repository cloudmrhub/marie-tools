function [Z] = Assembly_WSIE_block_par(WIRE,SCOIL,emc,inputs) 

    index_coil = SCOIL.index;
    etod_coil  = SCOIL.etod;
    node_coil  = SCOIL.node;
    elem_coil  = SCOIL.elem;
    first_p    = WIRE.F_point;
    second_p   = WIRE.S_point;
    third_p    = WIRE.T_point;
    
    [U,Z.V] = assembly_wire_surf_ns_aca_fun(first_p,second_p,third_p,index_coil,etod_coil,node_coil,elem_coil,inputs.Quad_order_wie,inputs.N_NS,emc.k0,inputs.tol_ACA);        
    Z.U = -(emc.eta0 / (4 * pi)) * U;

end