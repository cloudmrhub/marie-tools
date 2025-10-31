function I_S = surface_surface_kernels_Kop(I1_co,I2_co,I3_co,I4_co,W,A,B,C,D,Np,k0,r_m,r_n,m,dx,L)

    % 4D jacobean
    J = (dx/2)^4;

    % normal vector components
    nx = [-1 1  0 0  0 0];
    ny = [ 0 0 -1 1  0 0];
    nz = [ 0 0  0 0 -1 1];

    % n,n',l,l',pq
    I1 = zeros(6,6,L,L,3);
    % n,n',l,l',pq
    I2 = zeros(6,6,L,L,3);
    % n,n',l,l',pq
    I3 = zeros(6,6,L,L,3);
    % n,n',l,l',pq
    I4 = zeros(6,6,L,L,3);


    % calculate the adjacency type (ST,EA,VA) of squares of the square
    % of the faces of the voxel-voxel and the (4,6,7) points in proper
    % order for the calculation of the singular integrals
    [ adjacency_type , ordered_points ] = points_mapping(m,dx);

    for face = 1:6
        for face_p = 1:6

            n  = zeros(3,1);
            np = zeros(3,1);

            % calculate the normal vectors
            n (1)  = nx(face);
            np(1)  = nx(face_p);
            n (2)  = ny(face);
            np(2)  = ny(face_p);
            n (3)  = nz(face);
            np(3)  = nz(face_p);

            % calculate the centre of the squares
            r_ms = r_m + n *dx/2;
            r_ns = r_n + np*dx/2;

            % calculate the points for all 6 variables
            [X,Y,Z,Xp,Yp,Zp] = points_const_4D(A,B,C,D,r_ms,r_ns,dx,n,np);

            % get the points for the specific square-square
            R = ordered_points(:,:,face,face_p);


            for kernel = 1:4

                switch kernel 

                    case 1
                        %% surface - surface kernel 1
                        for l = 1:L
                            for lp = 1:L

                                switch adjacency_type(face,face_p)

                                    case 0
                                        % non-singular square-square 
                                        %evaluate kernel at the points of interest
                                        f1 = kernels_Kop(X,Y,Z,Xp,Yp,Zp,k0,dx,r_m,r_n,l,lp,n,1);

                                    case 1
                                        % self term
                                        fS1 = singular_ST_lin(l-1,lp-1,k0,dx,Np,5,n,np,R(:,1),R(:,2),R(:,3),R(:,4));

                                    case 2
                                        % edge adjacent term
                                        fS1 = singular_EA_lin(l-1,lp-1,k0,dx,Np,5,n,np,R(:,1),R(:,2),R(:,3),R(:,4),R(:,5),R(:,6));

                                    case 3
                                        % vertex adjacent term
                                        fS1 = singular_VA_lin(l-1,lp-1,k0,dx,Np,5,n,np,R(:,1),R(:,2),R(:,3),R(:,4),R(:,5),R(:,6),R(:,7));

                                end

                                 for i = 1:3

                                     % zero coefficient
                                     if I1_co(face,face_p,i) ~= 0
                                         % non-singular
                                         if adjacency_type(face,face_p) == 0
                                            I1(face,face_p,l,lp,i) = I1_co(face,face_p,i) * sum( W .*f1 ) * J;
                                         % singular
                                         else
                                            I1(face,face_p,l,lp,i) = I1_co(face,face_p,i) * fS1;
                                         end
                                     end

                                 end

                            end
                        end



                    case 2
                        %% surface - surfcae kernel 2
                        switch adjacency_type(face,face_p)

                                    case 0
                                        % non-singular square-square 
                                        %evaluate kernel at the points of interest
                                        f2 = kernels_Kop(X,Y,Z,Xp,Yp,Zp,k0,dx,r_m,r_n,0,0,n,2);

                                    case 1
                                        % self term
                                        fS2 = singular_ST_lin(0,0,k0,dx,Np,6,n,np,R(:,1),R(:,2),R(:,3),R(:,4));

                                    case 2
                                        % edge adjacent term
                                        fS2 = singular_EA_lin(0,0,k0,dx,Np,6,n,np,R(:,1),R(:,2),R(:,3),R(:,4),R(:,5),R(:,6));

                                    case 3
                                        % vertex adjacent term
                                        fS2 = singular_VA_lin(0,0,k0,dx,Np,6,n,np,R(:,1),R(:,2),R(:,3),R(:,4),R(:,5),R(:,6),R(:,7));
                        end

                        for lp = 1:L
                            for l = 1:L
                                for i = 1:3

                                    % zero coefficient
                                    if I2_co(face,face_p,i,l,lp) ~= 0
                                        % non-singular
                                        if adjacency_type(face,face_p) == 0 
                                            I2(face,face_p,l,lp,i) = I2_co(face,face_p,i,l,lp) * sum( W .*f2 ) * J;
                                        % singular    
                                        else
                                            I2(face,face_p,l,lp,i) = I2_co(face,face_p,i,l,lp) * fS2;
                                        end
                                    end

                                end
                            end

                        end

                    case 3
                        %% surface - surface kernel 3
                        for l = 1:L

                            switch adjacency_type(face,face_p)

                                    case 0
                                        % non-singular square-square 
                                        %evaluate kernel at the points of interest
                                        f3 = kernels_Kop(X,Y,Z,Xp,Yp,Zp,k0,dx,r_m,r_n,l,0,n,3);

                                    case 1
                                        % self term
                                        fS3 = singular_ST_lin(l-1,0,k0,dx,Np,7,n,np,R(:,1),R(:,2),R(:,3),R(:,4));

                                    case 2
                                        % edge adjacent term
                                        fS3 = singular_EA_lin(l-1,0,k0,dx,Np,7,n,np,R(:,1),R(:,2),R(:,3),R(:,4),R(:,5),R(:,6));

                                    case 3
                                        % vertex adjacent term
                                        fS3 = singular_VA_lin(l-1,0,k0,dx,Np,7,n,np,R(:,1),R(:,2),R(:,3),R(:,4),R(:,5),R(:,6),R(:,7));

                            end

                            for lp = 1:L
                                for i = 1:3

                                    % zero coefficient
                                    if I3_co(face,face_p,i,lp) ~= 0 
                                        % non-singular
                                        if adjacency_type(face,face_p) == 0 
                                            I3(face,face_p,l,lp,i) = I3_co(face,face_p,i,lp) * sum( W .*f3 ) * J;
                                        % singular    
                                        else
                                            I3(face,face_p,l,lp,i) = I3_co(face,face_p,i,lp) * fS3;
                                        end
                                    end

                                end
                            end

                        end

                    case 4 
                        %% surface - surface kernel 4

                        for lp = 1:L

                            switch adjacency_type(face,face_p)

                                        case 0
                                            % non-singular square-square 
                                            %evaluate kernel at the points of interest
                                            f4 = kernels_Kop(X,Y,Z,Xp,Yp,Zp,k0,dx,r_m,r_n,0,lp,n,4);

                                        case 1
                                            % self term
                                            fS4 = singular_ST_lin(0,lp-1,k0,dx,Np,8,n,np,R(:,1),R(:,2),R(:,3),R(:,4));

                                        case 2
                                            % edge adjacent term
                                            fS4 = singular_EA_lin(0,lp-1,k0,dx,Np,8,n,np,R(:,1),R(:,2),R(:,3),R(:,4),R(:,5),R(:,6));

                                        case 3
                                            % vertex adjacent term
                                            fS4 = singular_VA_lin(0,lp-1,k0,dx,Np,8,n,np,R(:,1),R(:,2),R(:,3),R(:,4),R(:,5),R(:,6),R(:,7));

                            end

                            for l = 1:L
                                 for i = 1:3

                                        % zero coefficient
                                        if I4_co(face,face_p,i,l) ~= 0
                                            % non-singular
                                            if adjacency_type(face,face_p) == 0
                                                I4(face,face_p,l,lp,i) = I4_co(face,face_p,i,l) * sum( W .*f4 ) * J;
                                            % singular
                                            else
                                                I4(face,face_p,l,lp,i) = I4_co(face,face_p,i,l) * fS4;
                                            end
                                        end

                                 end
                            end

                        end

                end
            end

        end
    end


    %% sum over 36 faces-faces
    I_SS = zeros(L,L,3);
    for l = 1:L
        for lp = 1:L
            for i = 1:3

                I_SS(l,lp,i) = sum(sum( I1(:,:,l,lp,i) + I2(:,:,l,lp,i) + I3(:,:,l,lp,i) + I4(:,:,l,lp,i) ));

            end
        end
    end

    % return only the unique integrals 
    I_S(1,:) = I_SS(1,1,:);

    if L == 4
        I_S(2,:) = I_SS(2,2,:);
        I_S(3,:) = I_SS(3,3,:);
        I_S(4,:) = I_SS(4,4,:);

        I_S(5,:) = I_SS(1,2,:);
        I_S(6,:) = I_SS(1,3,:);
        I_S(7,:) = I_SS(1,4,:);

        I_S(8,:)  = I_SS(2,3,:);
        I_S(9,:)  = I_SS(2,4,:);
        I_S(10,:) = I_SS(3,4,:);
    end

end