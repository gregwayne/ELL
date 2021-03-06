function [prototypes,bfns] = compute_granule_prototypes(K,N)

    % K number of prototypes
    % N number of principal components to use

    [celltypes,bfns]    = generate_bases('raw');

    % compute principal components
    [signals,PC,V]      = shlens_pca(bfns');
    wh_proj_sigs        = signals(1:N,:)./repmat(sqrt(V(1:N)),1,size(signals,2));

    % compute k-medoids
    cl_assns    = zeros(size(bfns,1),1);
    cl_meds     = wh_proj_sigs(:,1:K); % initialize with first K signals

    max_iter    = 50;

    for i=1:max_iter

        % assign to clusters
        for j=1:size(bfns,1)

            sig = wh_proj_sigs(:,j);

            min_diff = inf;
            for m=1:K

                cur_diff        = norm(sig-cl_meds(:,m));
                if cur_diff < min_diff
                    min_diff    = cur_diff;
                    cl_assns(j) = m;
                end

            end

        end
        
        % recompute medoids
        % swap in every other data point and take that as the cluster medoid 
        % if it has lowered cost
        for j=1:K

            midxs               = find(cl_assns==j);
            min_cost            = inf;
            min_m               = 0;
            for p=1:length(midxs)
                m=midxs(p);
                cst = medoid_cost(cl_assns,wh_proj_sigs(:,m),wh_proj_sigs,midxs);
                if cst < min_cost
                    min_cost    = cst;
                    min_m       = m;
                end

            end
            
            cl_meds(:,j) = wh_proj_sigs(:,min_m);
        end

    end
    
    prototypes = zeros(K,1);
    for i=1:K
       
        med = cl_meds(:,i);           
        for j=1:size(bfns,1)
           
            if norm(wh_proj_sigs(:,j)-med) == 0
                prototypes(i) = j;
            end
            
        end
        
    end    
    
end

function cst = medoid_cost(cl_assns,medoid,wh_proj_sigs,midxs)

    cst = 0;
    for p=1:length(midxs)
        i=midxs(p);

        cst = cst + norm(wh_proj_sigs(:,i)-medoid);
        
    end
        
end