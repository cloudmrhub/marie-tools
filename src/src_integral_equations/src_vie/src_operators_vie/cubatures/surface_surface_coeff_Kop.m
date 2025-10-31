function [I1_co,I2_co,I3_co,I4_co] = surface_surface_coeff_Kop(dx,ko,L)

    % normal vector components
    nx = [-1 1  0 0  0 0];
    ny = [ 0 0 -1 1  0 0];
    nz = [ 0 0  0 0 -1 1];

    % since G_K is anti-symmetric it suffices to calculate 3 unique interactions
    % p = [z x y]
    % q = [y z x]
      p = [3 1 2];
      q = [2 3 1];


    % n,n',pq
    I1_co = zeros(6,6,3);
    % n,n',pq,l,l'
    I2_co = zeros(6,6,3,L,L);
    % n,n',pq,l'
    I3_co = zeros(6,6,3,L);
    % n,n',pq,l
    I4_co = zeros(6,6,3,L);

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

            for coeff = 1:4

                switch coeff
                    case 1
                        %%
                        for i = 1:3
                            % calculate constant coefficient 1
                            I1_co(face,face_p,i) = coefficients_Kop(dx,ko,n,np,p(i),q(i),0,0,1);
                        end

                    case 2
                        %%
                        for i = 1:3
                            for l = 1:L
                                for lp = 1:L
                                    % calculate constant coefficient 2
                                    I2_co(face,face_p,i,l,lp) = coefficients_Kop(dx,ko,n,np,p(i),q(i),l,lp,2);
                                end
                            end
                        end

                    case 3
                        %%
                        for i = 1:3
                            for lp = 1:L
                                % calculate constant coefficient 3
                                I3_co(face,face_p,i,lp) = coefficients_Kop(dx,ko,n,np,p(i),q(i),0,lp,3);
                            end
                        end
                    case 4
                        %%
                        for i = 1:3
                            for l = 1:L
                                % calculate constant coefficient 4
                                I4_co(face,face_p,i,l) = coefficients_Kop(dx,ko,n,np,p(i),q(i),l,0,4);
                            end
                        end
                end

            end

        end
    end
  
 
end