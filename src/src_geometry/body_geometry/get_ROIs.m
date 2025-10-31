function [ROI] = get_ROIs(RHBM)

    x = RHBM.r(:,:,:,1);
    y = RHBM.r(:,:,:,2);
    z = RHBM.r(:,:,:,3);
    [n1,n2,n3] = size(x);
    body_mask = zeros(n1,n2,n3);
    body_mask(RHBM.idxS) = 1;
    
    rad = 0.04;
    zbot = 0.59;
    ztop = 0.78;
    x0 = -0.005;
    y0 = 0.02;
    ROI(1).id = find( sqrt((x - x0).^2 + (y - y0).^2) <= rad & z >= zbot & z <= ztop );
    ROI(2).id = find( RHBM.sigma_e>0.971 & RHBM.sigma_e<0.973 & z >= zbot);
    
    x1 = -0.04;
    x2 = -0.175;
    y1 = 0.035;
    y2 = 0.035;
    z1 = 0.62;
    z2 = 0.55;
    R = 0.035;
    vx = x2 - x1;  
    vy = y2 - y1;  
    vz = z2 - z1;
    vv = vx^2 + vy^2 + vz^2;
    ROI(3).id = find( ...
                    ( ( (x-x1)*vx + (y-y1)*vy + (z-z1)*vz ) ./ vv ) >= 0 & ... 
                    ( ( (x-x1)*vx + (y-y1)*vy + (z-z1)*vz ) ./ vv ) <= 1 & ...
                    ( (x-x1 - ( ( (x-x1)*vx + (y-y1)*vy + (z-z1)*vz )./vv )*vx ).^2 ... 
                    + (y-y1 - ( ( (x-x1)*vx + (y-y1)*vy + (z-z1)*vz )./vv )*vy ).^2 ...
                    + (z-z1 - ( ( (x-x1)*vx + (y-y1)*vy + (z-z1)*vz )./vv )*vz ).^2 ) <= R^2 );  
    
    x1 = 0.03;
    x2 = 0.155;
    y1 = 0.035;
    y2 = 0.045;
    z1 = 0.62;
    z2 = 0.55;
    R = 0.035;
    vx = x2 - x1;  
    vy = y2 - y1;  
    vz = z2 - z1;
    vv = vx^2 + vy^2 + vz^2;
    ROI(4).id = find( ...
                    ( ( (x-x1)*vx + (y-y1)*vy + (z-z1)*vz ) ./ vv ) >= 0 & ... 
                    ( ( (x-x1)*vx + (y-y1)*vy + (z-z1)*vz ) ./ vv ) <= 1 & ...
                    ( (x-x1 - ( ( (x-x1)*vx + (y-y1)*vy + (z-z1)*vz )./vv )*vx ).^2 ... 
                    + (y-y1 - ( ( (x-x1)*vx + (y-y1)*vy + (z-z1)*vz )./vv )*vy ).^2 ...
                    + (z-z1 - ( ( (x-x1)*vx + (y-y1)*vy + (z-z1)*vz )./vv )*vz ).^2 ) <= R^2 );      
    
    dist_to_air = bwdist(~body_mask);   
    core_mask   = dist_to_air > 3; 
    idx_core = find(core_mask);
    
    for k = 1:length(ROI)
        ROI(k).id = intersect(ROI(k).id, idx_core);  
    end

end

% all_mask = zeros(n1,n2,n3);
% all_mask(ROI(1).id) = 1;
% all_mask(ROI(2).id) = 1;
% all_mask(ROI(3).id) = 1;
% all_mask(ROI(4).id) = 1;
% cut_c = 30:42;
% 
% figure(1)
% fig_len = length(cut_c);
% c = 0;
% for i = 1:fig_len
%     c = c + 1;
%     subplot(fig_len,2,c)
%     imagesc(rot90(squeeze(RHBM.sigma_e(:,cut_c(i),:))))
%     colormap hot
%     colorbar
%     axis image
%     axis off
%     impixelinfo
%     clim([0 2.2])
% 
%     c = c + 1;
%     subplot(fig_len,2,c)
%     imagesc(rot90(squeeze(all_mask(:,cut_c(i),:).*RHBM.sigma_e(:,cut_c(i),:))))
%     colormap hot
%     colorbar
%     axis image
%     axis off
%     impixelinfo
%     clim([0 2.2])
% end