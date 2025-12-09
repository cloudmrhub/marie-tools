function y = prec_Block_Jacobi_LU_mult(Lprec, Pprec, Uprec, lims, x)

    coil = numel(lims)-1;
    y = zeros(lims(end),1,'like',x);
    s = lims(1:end-1)+1; 
    e = lims(2:end);

    for k = 1:coil
        Ik    = s(k):e(k);
        rhs   = x(Ik);
        rhsP  = rhs(Pprec(k).P);          
        tk    = Lprec(k).L * rhsP;         
        y(Ik) = Uprec(k).U * tk;           
    end
end
