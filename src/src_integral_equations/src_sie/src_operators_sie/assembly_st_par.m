function [ZR,ZR_DE_losses_matrix] = assembly_st_par(index,etod,node,elem,ZR,GL_order,Index_elem,emc) 

ko= emc.k0;
coil_losses = emc.rho_s;
eta = emc.eta0;
ZR_DE_losses_matrix = sparse(size(ZR,1),size(ZR,2));
%%%%%%%%%%%%%%% Gauss quadrature rule for singular triangles %%%%%%%%%%%%%%
Np_1D_WS_ST = GL_order.ST;
[w,z] = gauss_1d(Np_1D_WS_ST);
w = w';
z = z';
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Edge Adjacent Terms %%%%%%%%%%%%%%%%%%%%%%%%%%
ie_ST = Index_elem.ST(:,1); 
% je_NS = Index_elem.NS(:,2); 
n_ST_elem = length(ie_ST);

% semivectorized form
Z_ST_vector = zeros(9,n_ST_elem);
Z_ST_vector_losses = zeros(9,n_ST_elem);
AO_index_vector = zeros(9,n_ST_elem);
AS_index_vector = zeros(9,n_ST_elem);

parfor index_ST = 1 : n_ST_elem
    ie = ie_ST(index_ST);   % je = je_ST(index_ST);
    %
    ao = abs(etod(:,ie));
    % coordinates of nodes of the observation triangle
    node_test_1 = elem(1,ie);
    node_test_2 = elem(2,ie);
    node_test_3 = elem(3,ie);
    %
    ro_1 = node(:,node_test_1);
    ro_2 = node(:,node_test_2);
    ro_3 = node(:,node_test_3);
    %--------------------------------------------------------
    %                        DEMCEM                         %  
    [I_DE] = direct_ws_st_rwg(ro_1,ro_2,ro_3,ko,Np_1D_WS_ST,w,z);
    %--------------------------------------------------------
    as = ao;
    
    Z_ST_local = zeros(9,1);
    Z_ST_local_losses = zeros(9,1);
    AO_index_local = zeros(9,1);
    AS_index_local = zeros(9,1);
    localindex = 0;
            
    lo = [ro_2-ro_3,ro_3-ro_1,ro_1-ro_2];         % Losses
    aa = sqrt(lo(1,1)^2 + lo(2,1)^2 + lo(3,1)^2); % Losses
    bb = sqrt(lo(1,2)^2 + lo(2,2)^2 + lo(3,2)^2); % Losses
    cc = sqrt(lo(1,3)^2 + lo(2,3)^2 + lo(3,3)^2); % Losses
    s = (aa+bb+cc)/2;                             % Losses
    As = sqrt(s*(s-aa)*(s-bb)*(s-cc));            % Losses
    
    for i=1:3
        for j=1:3
            index_ao = index(ao(i));
            index_as = index(as(j));

            ii1=mod((i-1)+1,3)+1;                        % Losses
            ii2=mod((i-1)+2,3)+1;                        % Losses
            jj1=mod((j-1)+1,3)+1;                        % Losses
            jj2=mod((j-1)+2,3)+1;                        % Losses
            lii=sqrt(lo(1,i)^2 + lo(2,i)^2 + lo(3,i)^2); % Losses
            ljj=sqrt(lo(1,j)^2 + lo(2,j)^2 + lo(3,j)^2); % Losses
            
            if (index_ao && index_as)
                % sign of the dofs
                soi=sign(etod(i,ie));
                ssj=sign(etod(j,ie));
                %
                ZR_DE = soi*ssj*I_DE(i + 3*(j-1));
                % Losses -------------------------------------------------- 
                staticq = 0;
                if(i==j)
                    staticq = ( dot(lo(:,ii1),lo(:,ii1))/12 + dot(lo(:,ii2),lo(:,ii2))/12 - dot(lo(:,ii1),lo(:,ii2))/12);
                elseif(i==jj1)
                    staticq = (-dot(lo(:,jj1),lo(:,j))/12 + dot(lo(:,jj2),lo(:,j))/24 - dot(lo(:,jj2),lo(:,jj2))/24 + dot(lo(:,jj2),lo(:,jj1))/24);
                elseif(i==jj2)
                    staticq = (-dot(lo(:,jj2),lo(:,j))/12 + dot(lo(:,jj1),lo(:,j))/24 - dot(lo(:,jj1),lo(:,jj1))/24 + dot(lo(:,jj1),lo(:,jj2))/24);
                end
                ZR_DE_losses = 4*pi/eta*coil_losses*1/(2*As)*soi*ssj*lii*ljj*staticq;
                % ---------------------------------------------------------  
                
                localindex = localindex + 1;
                Z_ST_local(localindex,1) = ZR_DE+ZR_DE_losses;
                Z_ST_local_losses(localindex,1) = ZR_DE_losses;
                AO_index_local(localindex,1) = index_ao;
                AS_index_local(localindex,1) = index_as;
            end
            
        end
    end
    
    Z_ST_vector(:,index_ST) = Z_ST_local+Z_ST_local_losses;
    Z_ST_vector_losses(:,index_ST) = Z_ST_local_losses;
    AO_index_vector(:,index_ST) = AO_index_local;
    AS_index_vector(:,index_ST) = AS_index_local;
    
end


aoidx = nonzeros(AO_index_vector);
asidx = nonzeros(AS_index_vector);
zvals = nonzeros(Z_ST_vector);
zvals_losses = nonzeros(Z_ST_vector_losses);

for ii = 1:length(aoidx)
    ZR(aoidx(ii),asidx(ii)) = ZR(aoidx(ii),asidx(ii)) + zvals(ii);
end
for ii = 1:length(aoidx)
    ZR_DE_losses_matrix(aoidx(ii),asidx(ii)) = ZR_DE_losses_matrix(aoidx(ii),asidx(ii)) + zvals_losses(ii);
end
