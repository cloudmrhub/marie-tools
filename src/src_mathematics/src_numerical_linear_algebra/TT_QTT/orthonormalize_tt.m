function [tt_out,C4] = orthonormalize_tt(tt)

    n    = tt.n; 
    ps   = tt.ps; 
    core = tt.core; 
    r    = tt.r;
    
    g1      = core(ps(1):ps(2)-1);
    G1      = reshape(g1,n(1),r(2));
    [Q1,R1] = qr(G1,0);

    g2      = core(ps(2):ps(3)-1);
    G2      = reshape(g2,[r(2),n(2)*r(3)]);
    R1G2    = R1*G2;
    R1G2    = reshape(R1G2,[r(2)*n(2),r(3)]);
    [Q2,R2] = qr(R1G2,0);

    g3      = core(ps(3):ps(4)-1);
    G3      = reshape(g3,[r(3),n(3)*r(4)]);
    R2G3    = R2*G3;
    R2G3    = reshape(R2G3,[r(3)*n(3),r(4)]);
    [Q3,R3] = qr(R2G3,0);

    g4      = core(ps(4):ps(5)-1);
    G4      = reshape(g4,[r(4),n(4)]);
    R3G4    = R3*G4;
    R3G4    = reshape(R3G4,[r(4),n(4)]);
    [U,S,V] = svd(R3G4,'econ');


    core = zeros(length(ps(1):ps(4)-1),1);
    core(ps(1):ps(2)-1) = Q1(:);
    core(ps(2):ps(3)-1) = Q2(:);
    core(ps(3):ps(4)-1) = Q3(:);
    C4.U = U;
    C4.S = S;
    C4.V = V;

    tt_out.core = core;
    tt_out.ps = ps(1:4);
    tt_out.n = n(1:3);
    tt_out.r = r(1:4);

end
