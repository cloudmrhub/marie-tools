function [edge,etod,index,port,index_elem] = Mesh_PreProc(e,elem,port)
    % Pre Processing for GMSH discretization

    Nbd = size(e,2);    % Number of exterior (boundary) edges
    Ne  = size(elem,2); % Number of elements

    first_node  = [3 1 2]; % define edge vectors clockwise
    second_node = [2 3 1];

    % for edge and etod generation
    edge = zeros(2,3*Ne); % allocate maximun possible size
    etod = zeros(3,Ne); % allocate size

    % define array for boundary
    kn=zeros(3*Ne,1); % allocate maximun possible size

    % define arrays for adjacency
    eparent = zeros(3*Ne,1); % allocate maximum possible size for adjacency check

    NS = zeros(floor(Ne*Ne/2),2); nscount = 0;
    VA = zeros(100*Ne,2); vacount = 0;
    EA = zeros(3*Ne,2); eacount = 0;

    % loop on elements
    Nd = 0;
    for ii=1:Ne
        count = Nd;
        flagpad = zeros(ii-1,1);
        flagead = zeros(ii-1,1);

        for in=1:3
            node1 = elem(first_node(in),ii);        
            node2 = elem(second_node(in),ii);
            flage  = 0;
            R = elem(1,1:ii-1) - node1; flagpad(R == 0) = 1;
            R = elem(2,1:ii-1) - node1; flagpad(R == 0) = 1;
            R = elem(3,1:ii-1) - node1; flagpad(R == 0) = 1;

            for jj=1:Nd
                if (node1 == edge(1,jj)) && (node2 == edge(2,jj))
                    etod(in,ii) = jj; % positive match found
                    flage = 1; % edge adjacency
                    flagead(eparent(jj)) = 2;
                    % update kn if kn(jj) was exterior edge
                    if 0 == kn(jj)
                        kn(jj) = -1; % set boundary flag as internal node
                    end
                end
                if (node2 == edge(1,jj)) &&  (node1 == edge(2,jj))
                    etod(in,ii) = -jj; % negative match found
                    flage = 1; % edge adjacency
                    flagead(eparent(jj)) = 2;
                    % update kn if kn(jj) was exterior edge
                    if 0 == kn(jj)
                        kn(jj) = -1; % set boundary flag as internal node
                    end
                end
            end

            if (flage == 0) % new edge
                count = count+1;
                edge(1,count) = node1;
                edge(2,count) = node2;
                etod(in,ii) = count;
                eparent(count) = ii; % store the element to which the edge belongs
                for kk=1:Nbd % assign boundary
                    bn1=e(1,kk);
                    bn2=e(2,kk);
                    if( (bn1 == node1) && (bn2 == node2) )
                        kn(count)=e(3,kk); % zero if external, positive integer if port
                    end
                    if( (bn2 == node1) && (bn1 == node2))
                        kn(count)=e(3,kk); % zero if external, positive integer if port
                    end
                end
            end
        end
        %check adjacency
        for kk = 1:length(flagead)
            adval = max(flagead(kk),flagpad(kk));
            switch adval
                case 1
                    vacount = vacount+1;
                    VA(vacount,1) = ii;
                    VA(vacount,2) = kk;
                case 2
                    eacount = eacount+1;
                    EA(eacount,1) = ii;
                    EA(eacount,2) = kk;
                case 0
                    nscount = nscount+1;
                    NS(nscount,1) = ii;
                    NS(nscount,2) = kk;
            end
        end
        Nd = count;
    end
    edge = edge(:,1:Nd);
    kn = kn(1:Nd);
    % remove replicas in adjacency
    ST = [1:Ne; 1:Ne].';
    EA = EA(1:eacount,:);
    VA = VA(1:vacount,:);
    NS = NS(1:nscount,:);
	EA = sort(EA,2,"descend");
	VA = sort(VA,2,"descend");
	EA = sortrows(EA);
	VA = sortrows(VA);

    % map NS matrix to cell arrays
    NS = sort(NS,2);
    [~, NS_idx] = sort(NS(:,1));
    NS = NS(NS_idx,:);

    % Indexing & Number of unknown dofs (Nff)
    index=zeros(Nd,1);
    etype = sort(unique(kn)); 
    idx = etype > 0;
    etype = etype(idx);

    tE_list_mask = strcmp({port.type}, 'port');
    tL_list_mask = strcmp({port.type}, 'element');
    tE_list = etype(tE_list_mask);
    tL_list = etype(tL_list_mask);

    % loop over excitation ports first
    term_count = 0;
    for i = 1:length(tE_list)
        idx = find(kn == tE_list(i));    
        p_num = etype == tE_list(i);    
        dofnum = term_count+1:term_count+length(idx); 
        index(idx) = dofnum; 
        port(p_num).t = dofnum;      
        term_count = term_count+length(idx); 
    end

    % loop over load ports second
    for i = 1:length(tL_list)
        idx = find(kn == tL_list(i));   
        p_num = etype == tL_list(i);    
        dofnum = term_count+1:term_count+length(idx);     
        index(idx) = dofnum;    
        port(p_num).t = dofnum; 
        term_count = term_count+length(idx); 
    end

    % add rest of the edges
    idx = find(kn == -1);
    dofnum = term_count+1:term_count+length(idx); 
    index(idx) = dofnum;
    index_elem.ST = ST;
    index_elem.VA = VA;
    index_elem.EA = EA;
    index_elem.NS = NS;

end