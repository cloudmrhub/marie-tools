function [I1_co,I2_co,I3_co,I4_co] = surface_surface_coeff_Nop(dx,L)

    %% inputs and outputs 

    % dx    : resolution

    % I1_co : (6x6x6) double array with the constant coefficients that we can take out of the surface-surface integral for the 1st kernel
    % I2_co : (6x6x6x4)                                        -- // --                                                        2nd kernel
    % I3_co : (6x6x6x4)                                        -- // --                                                        3rd kernel
    % I4_co : (6x6x6x4x4)                                      -- // --                                                        4th kernel


    %% calculation of all the coefficients

    % normal vector components
    nx = [-1 1  0 0  0 0];
    ny = [ 0 0 -1 1  0 0];
    nz = [ 0 0  0 0 -1 1];

    % since G_N is symmetric it suffices to calculate 6 interactions
    % p = [x x x y y z]
    % q = [x y z y z z]
      p = [1 1 1 2 2 3];
      q = [1 2 3 2 3 3];

    % n,n',pq
    I1_co = zeros(6,6,6);
    % n,n',pq,l
    I2_co = zeros(6,6,6,L);
    % n,n',pq,l'
    I3_co = zeros(6,6,6,L);
    % n,n',pq,l,l'
    I4_co = zeros(6,6,6,L,L);

    % loop over the 36 squares-squares 
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

            % loop over the 4 surface-surface kernels
            for coeff = 1:4

                switch coeff
                    case 1
                        %%
                        for i = 1:6
                            % calculate constant coefficient for kernel 1
                            I1_co(face,face_p,i) = coefficients_Nop(dx,n,np,p(i),q(i),0,0,1);
                        end

                    case 2
                        %%
                        for i = 1:6
                            for l = 1:L
                                % calculate constant coefficient for kernel 2
                                I2_co(face,face_p,i,l) = coefficients_Nop(dx,n,np,p(i),q(i),l,0,2);
                            end
                        end

                    case 3
                        %%
                        for i = 1:6
                            for lp = 1:L
                                % calculate constant coefficient for kernel 3
                                I3_co(face,face_p,i,lp) = coefficients_Nop(dx,n,np,p(i),q(i),0,lp,3);
                            end
                        end
                    case 4
                        %%
                        for i = 1:6
                            for l = 1:L
                                for lp = 1:L
                                    % calculate constant coefficient for kernel 4
                                    I4_co(face,face_p,i,l,lp) = coefficients_Nop(dx,n,np,p(i),q(i),l,lp,4);
                                end
                            end
                        end
                end

            end

        end
    end
  
 
end