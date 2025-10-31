function [y]=dmrg_cross_cpu(d,n,fun,eps)

    kickrank=2;
    nswp=5;
    y=[];
    maxr = 500; 

    if ( numel(n) == 1 )
       n=n*ones(d,1);
    end
    sz=n;
    if (isempty(y) )
        y=tt_rand(sz,d,2); 
    end
    elem = fun; 
    y=round(y,0); 
    ry=y.r;
    [y,rm]=qr(y,'rl');
    y=rm*y;

    swp=1;
    rmat=cell(d+1,1); 
    rmat{d+1}=1;
    rmat{1}=1;
    index_array{d+1}=zeros(0,ry(d+1)); 
    index_array{1}=zeros(ry(1),0);
    r1=1;
    
    for i=d:-1:2
        cr=y{i}; cr=reshape(cr,[ry(i)*n(i),ry(i+1)]);
        cr = cr*r1; cr=reshape(cr,[ry(i),n(i)*ry(i+1)]); cr=cr.';
        [cr,rm]=qr(cr,0);
        [ind]=maxvol2(cr); 
        ind_old=index_array{i+1};
        rnew=min(n(i)*ry(i+1),ry(i));
        ind_new=zeros(d-i+1,rnew);
        for s=1:rnew
           f_in=ind(s);
           w1=tt_ind2sub([ry(i+1),n(i)],f_in);
           rs=w1(1); js=w1(2);
           ind_new(:,s)=[js,ind_old(:,rs)'];
        end
        index_array{i}=ind_new;
        r1=cr(ind,:);
        cr=cr/r1; 
        r1=r1*rm;
        r1=r1.';

        cr=cr.'; 
        y{i}=reshape(cr,[ry(i),n(i),ry(i+1)]);
        cr=reshape(cr,[ry(i)*n(i),ry(i+1)]);
        cr=cr*rmat{i+1}; cr=reshape(cr,[ry(i),n(i)*ry(i+1)]);
        cr=cr.'; 
        [~,rm]=qr(cr,0);
        rmat{i}=rm; 
    end
    
    cr=y{1}; cr=reshape(cr,[ry(1)*n(1),ry(2)]);
    y{1}=reshape(cr*r1,[ry(1),n(1),ry(2)]); 
    not_converged = true;
    dir = 1;
    i=1;
    er_max=0;
    
    while ( swp < nswp && not_converged )
        cr1=y{i}; cr2=y{i+1};
        ind1=index_array{i};
        ind2=index_array{i+2};

        big_index = [ ...
          ind1(repmat((1:ry(i))', n(i)*n(i+1)*ry(i+2), 1),:), ...
          kron(repmat((1:n(i))', n(i+1)*ry(i+2), 1), ones(ry(i),1)), ...
          kron(repmat((1:n(i+1))', ry(i+2), 1), ones(ry(i)*n(i),1)), ...          
          ind2(:, kron((1:ry(i+2))', ones(ry(i)*n(i)*n(i+1),1)))' ...      
          ]; 

        score=elem(big_index);     

        score=reshape(score,[ry(i),n(i)*n(i+1)*ry(i+2)]);
        score=rmat{i}*score;
        ry(i)=size(score,1);
        score=reshape(score,[ry(i)*n(i)*n(i+1),ry(i+2)]);
        score=score*rmat{i+2}; 
        ry(i+2)=size(score,2);

        score=reshape(score,[ry(i)*n(i),n(i+1)*ry(i+2)]);

%         [u,s,v]=svds(score,maxr);
        [u,s,v]=svd(score,"econ");
        s=diag(s); 
        r=my_chop2(s,norm(s)*eps/sqrt(d-1)); 
        u=u(:,1:r); v=v(:,1:r); s=diag(s(1:r));

        if ( dir == 1 ) 
            v = v * s'; 
            ur=randn(size(u,1),kickrank);
            u=reort(u,ur);
            radd=size(u,2)-r;
            if ( radd > 0 )
                vr=zeros(size(v,1),radd);
                v=[v,vr];
            end
            r=r+radd;
        else
             u = u * s; 
             vr=randn(size(v,1),kickrank);
             v=reort(v,vr);
             radd=size(v,2)-r;
             if ( radd > 0 )
                 ur=zeros(size(u,1),radd);
                 u=[u,ur];
             end
             r=r+radd;
        end

        v=v';

        appr=reshape(cr1,[numel(cr1)/ry(i+1),ry(i+1)])*reshape(cr2,[ry(i+1),numel(cr2)/ry(i+1)]);
        appr=reshape(appr,[ry(i),n(i)*n(i+1)*ry(i+2)]);
        appr=rmat{i}*appr;
        appr=reshape(appr,[ry(i)*n(i)*n(i+1),ry(i+2)]);
        appr=appr*rmat{i+2}; 
        er_loc=norm(score(:)-appr(:))/norm(score(:));
        er_max=max(er_max,er_loc);
        ry(i+1)=r;

        u = reshape(u,[ry(i),n(i)*r]);
        u = rmat{i}\u; 
        v=reshape(v,[r*n(i+1),ry(i+2)]); 
        u=reshape(u,[ry(i)*n(i),ry(i+1)]);
        v=v/rmat{i+2}; v=reshape(v,[r,n(i+1)*ry(i+2)]);
        if ( dir == 1 ) 
            [u,rm]=qr(u,0); 
            ind=maxvol2(u); 
            r1=u(ind,:); 
            u=u/r1; y{i}=reshape(u,[ry(i),n(i),ry(i+1)]);
            r1=r1*rm; 
            v=r1*v; y{i+1}=reshape(v,[ry(i+1),n(i+1),ry(i+2)]);
            u1=reshape(u,[ry(i),n(i)*ry(i+1)]);
            u1=rmat{i}*u1;
            u1=reshape(u1,[ry(i)*n(i),ry(i+1)]);
            [~,rm]=qr(u1,0);
            rmat{i+1}=rm;
            ind_old=index_array{i};
            ind_new=zeros(ry(i+1),i);
            for s=1:ry(i+1)
                f_in=ind(s);
                w1=tt_ind2sub([ry(i),n(i)],f_in);
                rs=w1(1); js=w1(2);
                ind_new(s,:)=[ind_old(rs,:),js];
            end
            index_array{i+1}=ind_new; 
            if ( i == d - 1 ) 
                dir = -dir;
            else
                i=i+1;
            end
        else
            v=v.'; 
            [v,rm]=qr(v,0);
            ind=maxvol2(v);
            r1=v(ind,:);
            v=v/r1; v2=reshape(v,[n(i+1),ry(i+2),ry(i+1)]); y{i+1}=permute(v2,[3,1,2]);
            r1=r1*rm; r1=r1.';
            u=u*r1; y{i}=reshape(u,[ry(i),n(i),ry(i+1)]);
            v=v.'; 
            v=reshape(v,[ry(i+1)*n(i+1),ry(i+2)]);
            v=v*rmat{i+2};
            v=reshape(v,[ry(i+1),n(i+1)*ry(i+2)]); v=v.';
            [~,rm]=qr(v,0);
            rmat{i+1}=rm;
            ind_old=index_array{i+2};
            ind_new=zeros(d-i,ry(i+1));
            for s=1:ry(i+1)
                f_in=ind(s);
                w1=tt_ind2sub([n(i+1),ry(i+2)],f_in);
                rs=w1(2); js=w1(1);
                ind_new(:,s)=[js,ind_old(:,rs)'];
            end
            index_array{i+1}=ind_new;
            if ( i == 1 ) 
                dir=-dir;
                swp = swp + 1;
                if ( er_max < eps ) 
                    not_converged=false;
                    fprintf('\tConvergence error: %4.6f \n',er_max);
                else
                    fprintf('\tConvergence error: %4.6f \n',er_max);
                    er_max=0;
                end
            else
                i=i-1;
            end
        end
    end
    return

end
