function visualize_ICP(target_point,az,el,SK,VK,wUISNR,cUISNR,Substrate,RHBM)
    
    set(groot,'defaultAxesTickLabelInterpreter','latex');  
    fs = 20;
    first_node  = [3 1 2];
    second_node = [2 3 1];

    xd = RHBM.r(:,:,:,1);
    yd = RHBM.r(:,:,:,2);
    zd = RHBM.r(:,:,:,3);
    idxS = RHBM.idxS;
    x2 = xd(idxS);
    y2 = yd(idxS);
    z2 = zd(idxS);
    P = [x2 y2 z2];
    k = boundary(P,1);

    % Compute Euclidean distance to the target point
    leng = size(target_point,2);
    j_cart_p_0 = zeros(size(Substrate.elem,2),3,leng);
    j_cart_p_pi2 = zeros(size(Substrate.elem,2),3,leng);

    for ii = 1:leng
        dist = sqrt((xd - target_point(1,ii)).^2 + ...
                    (yd - target_point(2,ii)).^2 + ...
                    (zd - target_point(3,ii)).^2);


        % Find index of minimum distance
        [~, id_voxel1] = min(dist(:));
        id_voxel = find(idxS==id_voxel1);

        j_cart = zeros(size(Substrate.elem,2),3);
        sh = (unique(floor(logspace(0,log10(size(SK,1)),50))))';
        modes_basis = size(cUISNR,2);
        for i = 1:size(cUISNR,2)
            if abs(cUISNR(id_voxel,i))/abs(cUISNR(id_voxel,end)) < 0.95
                continue;
            else
                modes_basis = i;
                break;
            end
        end
        modes_basis = sh(modes_basis);
        w = wUISNR(:,id_voxel);
        Jc = VK(:,1:modes_basis)/SK(1:modes_basis,1:modes_basis)*w(1:modes_basis,:);
        for jj = 1:size(Substrate.elem,2)
            rp  = Substrate.Ct(:,jj);
            r1  = Substrate.node(:,Substrate.elem(1,jj));
            r2  = Substrate.node(:,Substrate.elem(2,jj));
            r3  = Substrate.node(:,Substrate.elem(3,jj));
            Ae  = Triangle_area(r1,r2,r3);
            Ae1 = Triangle_area(rp,r2,r3);
            Ae2 = Triangle_area(r1,rp,r3);
            Ae3 = Triangle_area(r1,r2,rp);
            zs  = [Ae1 Ae2 Ae3]./Ae;
            ed  = abs(Substrate.etod(:,jj));
            lo  = [r2(1)-r3(1) r3(1)-r1(1) r1(1)-r2(1); r2(2)-r3(2) r3(2)-r1(2) r1(2)-r2(2); r2(3)-r3(3) r3(3)-r1(3) r1(3)-r2(3)];
            for id = 1:3
                if Substrate.index(ed(id))~=0
                    i_1 = second_node(id);
                    i_2 = first_node(id);
                    L_i = sqrt(lo(1,id)^2 + lo(2,id)^2 + lo(3,id)^2);
                    j_cart(jj,:) = j_cart(jj,:) + (L_i/(2*Ae)) * (zs(i_2)*lo(:,i_1)' - zs(i_1)*lo(:,i_2)') * sign(Substrate.etod(id,jj)) * Jc(Substrate.index(ed(id)),1);
                end
            end
        end
    
        wt = 0;
        jx = real(j_cart(:,1)*exp(1i*wt));
        jy = real(j_cart(:,2)*exp(1i*wt));
        jz = real(j_cart(:,3)*exp(1i*wt));
        j_cart_p_0(:,:,ii) = [jx jy jz];
        wt = pi/2;
        jx = real(j_cart(:,1)*exp(1i*wt));
        jy = real(j_cart(:,2)*exp(1i*wt));
        jz = real(j_cart(:,3)*exp(1i*wt));
        j_cart_p_pi2(:,:,ii) = [jx jy jz];
    end


    x_body = P(:,1)*100;
    y_body = P(:,2)*100;
    z_body = P(:,3)*100;
    x_substrate = Substrate.Ct(1,:)'*100;
    y_substrate = Substrate.Ct(2,:)'*100;
    z_substrate = Substrate.Ct(3,:)'*100;
    plot_nodes = Substrate.node'*100;
    TR = triangulation(Substrate.elem(1:3,:)',plot_nodes);

    figure

    [n1,n2,n3] = size(xd);
    [v1,v2,v3] = ind2sub([n1 n2 n3], id_voxel1);

    for ii = 1:leng
        for i = 1:2
            subplot(leng,2,i+2*(ii-1))
            trisurf(k,x_body,y_body,z_body,'LineStyle','none','FaceColor',[224 177 164]/255,'FaceAlpha',0.5); hold on;
            scatter3(xd(v1,v2,v3)*100,yd(v1,v2,v3)*100,zd(v1,v2,v3)*100,100,[255 44 44]/255,'s','filled'); hold on;
            trisurf(TR,'FaceColor',[255 238 117]/255,'FaceAlpha',0.3,'LineStyle','none'); 
            if i == 1
                quiver3(x_substrate,y_substrate,z_substrate,j_cart_p_0(:,1,ii),j_cart_p_0(:,2,ii),j_cart_p_0(:,3,ii),5.5,'k'); hold off;
                if ii == 1
                    title('$\omega t =0$','Interpreter','latex','Fontsize',fs);
                end
            elseif i == 2
                quiver3(x_substrate,y_substrate,z_substrate,j_cart_p_pi2(:,1,ii),j_cart_p_pi2(:,2,ii),j_cart_p_pi2(:,3,ii),5.5,'k'); hold off;
                if ii == 1
                    title('$\omega t = \pi/2$','Interpreter','latex','Fontsize',fs);
                end
            end
            axis equal;
            axis off;
            grid off;
            ax = gca;
            ax.FontSize = fs;
	        view(az(ii),el(ii));
	        camlight
	        lighting gouraud;
        end
    end
end