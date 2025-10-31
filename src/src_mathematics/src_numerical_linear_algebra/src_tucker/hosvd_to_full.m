function[A] = hosvd_to_full(G,U1,U2,U3)
            
    A = nmp(nmp(nmp(G,U1,1),U2,2),U3,3);

end