function [ginv_ax,ginv_cr,ginv_sg] = em_g_factor(phimat,b1minus,mask,cut_ax,cut_cr,cut_sg,Rp,Rf)

    phimat(phimat == diag(phimat)) = real(diag(phimat));
    L = chol(phimat,'lower');

    b1_ax               = permute(squeeze(b1minus(:,:,cut_ax,:)),[3 2 1]);
    mask_ax             = permute(squeeze(mask(:,:,cut_ax)),[2 1]);
    [nc,np_ax,nf_ax]    = size(b1_ax);
    b1_ax               = reshape(b1_ax,[nc np_ax*nf_ax]);
    b1map_ax            = L\b1_ax(:,:);
    b1map_ax            = reshape(b1map_ax,nc,np_ax,nf_ax);
    
    b1_cr               = permute(squeeze(b1minus(:,cut_cr,:,:)),[3 2 1]);
    mask_cr             = permute(squeeze(mask(:,cut_cr,:)),[2 1]);
    [~,np_cr,nf_cr] = size(b1_cr);
    b1_cr               = reshape(b1_cr,[nc np_cr*nf_cr]);
    b1map_cr            = L\b1_cr(:,:);
    b1map_cr            = reshape(b1map_cr,nc,np_cr,nf_cr);
    
    b1_sg               = permute(squeeze(b1minus(cut_sg,:,:,:)),[3 2 1]);
    mask_sg             = permute(squeeze(mask(cut_sg,:,:)),[2 1]);
    [~,np_sg,nf_sg] = size(b1_sg);
    b1_sg               = reshape(b1_sg,[nc np_sg*nf_sg]);
    b1map_sg            = L\b1_sg(:,:);
    b1map_sg            = reshape(b1map_sg,nc,np_sg,nf_sg);

    g_ax = zeros(np_ax,nf_ax);
    g_cr = zeros(np_cr,nf_cr);
    g_sg = zeros(np_sg,nf_sg);

    for x=1:floor(nf_ax./Rf)
        for y=1:floor(np_ax./Rp)
            s_temp=squeeze(b1map_ax(:,y:floor(np_ax./Rp):np_ax,x:floor(nf_ax./Rf):nf_ax));  
            s = reshape(s_temp,[nc size(s_temp,2)*size(s_temp,3)]);
            g_ax(y:floor(np_ax./Rp):np_ax,x:floor(nf_ax./Rf):nf_ax) = ...
                reshape(sqrt(abs(diag(pinv(s'*s)).*diag(s'*s))),[size(s_temp,2) size(s_temp,3)]);   
        end
    end

    for x=1:floor(nf_cr./Rf)
        for y=1:floor(np_cr./Rp)
            s_temp=squeeze(b1map_cr(:,y:floor(np_cr./Rp):np_cr,x:floor(nf_cr./Rf):nf_cr));  
            s = reshape(s_temp,[nc size(s_temp,2)*size(s_temp,3)]);
            g_cr(y:floor(np_cr./Rp):np_cr,x:floor(nf_cr./Rf):nf_cr) = ...
                reshape(sqrt(abs(diag(pinv(s'*s)).*diag(s'*s))),[size(s_temp,2) size(s_temp,3)]);   
        end
    end

    for x=1:floor(nf_sg./Rf)
        for y=1:floor(np_sg./Rp)
            s_temp=squeeze(b1map_sg(:,y:floor(np_sg./Rp):np_sg,x:floor(nf_sg./Rf):nf_sg));  
            s = reshape(s_temp,[nc size(s_temp,2)*size(s_temp,3)]);
            g_sg(y:floor(np_sg./Rp):np_sg,x:floor(nf_sg./Rf):nf_sg) = ...
                reshape(sqrt(abs(diag(pinv(s'*s)).*diag(s'*s))),[size(s_temp,2) size(s_temp,3)]);   
        end
    end

    g_ax = g_ax.*mask_ax;
    ginv_ax = (1./g_ax).*mask_ax;

    g_cr = g_cr.*mask_cr;
    ginv_cr = (1./g_cr).*mask_cr;

    g_sg = g_sg.*mask_sg;
    ginv_sg = (1./g_sg).*mask_sg;

    ginv_ax = ginv_ax.';
    ginv_cr = ginv_cr.';
    ginv_sg = ginv_sg.';

end

