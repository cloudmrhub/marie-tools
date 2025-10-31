function [W,X,Y,Z,Xp,Yp,Zp] = weights_points(N,dim)

    [ w1d , x1d ] = gauss_1d(N);

    if dim == 1

        W = w1d;
        X = x1d;

    elseif dim == 2

        w2d = w1d*transpose(w1d);
        W = w2d(:);
        [X,Y] = meshgrid(x1d);
        X = X(:);
        Y = Y(:);

    elseif dim == 3

        w2d = w1d*transpose(w1d);
        w3d = w2d(:)*transpose(w1d);
        W = w3d(:);
        [X,Y,Z] = meshgrid(x1d);
        X = X(:);
        Y = Y(:);
        Z = Z(:);

    elseif dim == 4

        w2d = w1d*transpose(w1d);
        w3d = w2d(:)*transpose(w1d);
        w4d = w3d(:)*transpose(w1d);
        W = w4d(:);
        [X,Y,Z,Xp] = ndgrid(x1d);
        X = X(:);
        Y = Y(:);
        Z = Z(:);
        Xp = Xp(:);  

    elseif dim == 5

        w2d = w1d*transpose(w1d);
        w3d = w2d(:)*transpose(w1d);
        w4d = w3d(:)*transpose(w1d);
        w5d = w4d(:)*transpose(w1d);
        W = w5d(:);
        [X,Y,Z,Xp,Yp] = ndgrid(x1d);
        X = X(:);
        Y = Y(:);
        Z = Z(:);
        Xp = Xp(:);
        Yp = Yp(:);

    elseif dim == 6

        w2d = w1d*transpose(w1d);
        w3d = w2d(:)*transpose(w1d);
        w4d = w3d(:)*transpose(w1d);
        w5d = w4d(:)*transpose(w1d);
        w6d = w5d(:)*transpose(w1d);
        W = w6d(:);
        [X,Y,Z,Xp,Yp,Zp] = ndgrid(x1d);
        X = X(:);
        Y = Y(:);
        Z = Z(:);
        Xp = Xp(:);
        Yp = Yp(:);
        Zp = Zp(:);

    end    
    
end