function  I = VV_Kop(W,X,Y,Z,Xp,Yp,Zp,k0,r_m,r_n,dx,pwx)

    J = (dx/2)^6;

    Xg = r_m(1) + dx/2*X;
    Yg = r_m(2) + dx/2*Y;
    Zg = r_m(3) + dx/2*Z;

    Xpg = r_n(1) + dx/2*Xp;
    Ypg = r_n(2) + dx/2*Yp;
    Zpg = r_n(3) + dx/2*Zp;

    GK = G_K(Xg,Yg,Zg,Xpg,Ypg,Zpg,k0);

    Xt  = dx/2*X;
    Yt  = dx/2*Y;
    Zt  = dx/2*Z;

    Xpb  = dx/2*Xp;
    Ypb  = dx/2*Yp;
    Zpb  = dx/2*Zp;

    I = zeros(pwx,3);

    ind_l  = [1 2 3 4 1 1 1 2 2 3];
    ind_lp = [1 2 3 4 2 3 4 3 4 4];

     for llp = 1:pwx

         l = ind_l(llp);
         lp = ind_lp(llp);

         switch l
             case 1
                 T = 1;
             case 2
                 T = Xt/dx;
             case 3
                 T = Yt/dx;
             case 4
                 T = Zt/dx;
         end

         switch lp
             case 1
                 B = 1;
             case 2
                 B = Xpb/dx;
             case 3
                 B = Ypb/dx;
             case 4
                 B = Zpb/dx;
         end

         TB = T.*B;

         for i = 1:3

             kernel = TB .* GK(:,i);         
             I(llp,i) = sum ( W .* kernel ) * J;

         end

     end
 
end