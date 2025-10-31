function f = G_N(x,y,z,xp,yp,zp,k0)

    R = sqrt((x-xp).^2+(y-yp).^2+(z-zp).^2);
    Gexp = exp(-1i*k0*R)/4/pi;

    Gxx = Gexp .* ( 3*(x-xp).^2./R.^5 + 3i*k0*(x-xp).^2./R.^4 - (1+k0^2*(x-xp).^2)./R.^3 - 1i*k0./R.^2 );
    Gyy = Gexp .* ( 3*(y-yp).^2./R.^5 + 3i*k0*(y-yp).^2./R.^4 - (1+k0^2*(y-yp).^2)./R.^3 - 1i*k0./R.^2 );
    Gzz = Gexp .* ( 3*(z-zp).^2./R.^5 + 3i*k0*(z-zp).^2./R.^4 - (1+k0^2*(z-zp).^2)./R.^3 - 1i*k0./R.^2 );

    Gxy = Gexp .* ( 3*(x-xp).*(y-yp)./R.^5 + 3i*k0*(x-xp).*(y-yp)./R.^4 -k0^2*(x-xp).*(y-yp)./R.^3 );
    Gxz = Gexp .* ( 3*(x-xp).*(z-zp)./R.^5 + 3i*k0*(x-xp).*(z-zp)./R.^4 -k0^2*(x-xp).*(z-zp)./R.^3 );
    Gyz = Gexp .* ( 3*(y-yp).*(z-zp)./R.^5 + 3i*k0*(y-yp).*(z-zp)./R.^4 -k0^2*(y-yp).*(z-zp)./R.^3 );

    f(:,1) = -Gyy - Gzz;
    f(:,2) =  Gxy;
    f(:,3) =  Gxz;
    f(:,4) = -Gxx - Gzz;
    f(:,5) =  Gyz;
    f(:,6) = -Gxx - Gyy;

end

