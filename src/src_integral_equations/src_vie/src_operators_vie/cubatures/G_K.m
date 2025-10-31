function f = G_K(x,y,z,xp,yp,zp,k0)

    R = sqrt((x-xp).^2+(y-yp).^2+(z-zp).^2);
    Gexp = exp(-1i*k0*R)/4/pi;

    Gx = Gexp .* (-(x-xp)) .* (1./R.^3 + 1i*k0./R.^2);
    Gy = Gexp .* (-(y-yp)) .* (1./R.^3 + 1i*k0./R.^2);
    Gz = Gexp .* (-(z-zp)) .* (1./R.^3 + 1i*k0./R.^2);

    f(:,1) = Gx;
    f(:,2) = Gy;
    f(:,3) = Gz;

end

