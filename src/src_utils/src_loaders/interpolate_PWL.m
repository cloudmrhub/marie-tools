function [A_f] = interpolate_PWL(A,stepi)
    
    [L,M,N,~] = size(A);
    
    L_f = zeros(L,stepi);
    M_f = zeros(M,stepi);
    N_f = zeros(N,stepi);
    
    for i = 1:stepi
        L_f(:,i) = i:stepi:stepi*L-(stepi-i);
        M_f(:,i) = i:stepi:stepi*M-(stepi-i);
        N_f(:,i) = i:stepi:stepi*N-(stepi-i);
    end
    
    A_f = zeros(stepi*L,stepi*M,stepi*N,3);
    
    for l = 1:stepi
        for m = 1:stepi
            for n = 1:stepi
                
                l_step = -(stepi-(2*l-1))/(2*stepi);
                m_step = -(stepi-(2*m-1))/(2*stepi);
                n_step = -(stepi-(2*n-1))/(2*stepi);
                
                for i = 1:3   
                    
                    A_f (L_f(:,l),M_f(:,m),N_f(:,n),i) = A(:,:,:,(i-1)*4+1) +  l_step * A(:,:,:,(i-1)*4+2) + m_step *A(:,:,:,(i-1)*4+3) + n_step *A(:,:,:,(i-1)*4+4);
                    
                end
                
            end
        end
    end
    
end