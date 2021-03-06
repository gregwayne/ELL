function prototypes = compute_granule_prototypes()

    % todo: get the units on the protoypes right
    K = 4; % number of prototypes
    N = 5; % number of principal components to use

    [celltypes,bfns]    = generate_bases('raw');

    % compute principal components
    [signals,PC,V]      = shlens_pca(bfns');
    wh_proj_sigs        = signals(1:N,:)./repmat(sqrt(V(1:N)),1,size(signals,2));

    scatter3(wh_proj_sigs(:,1),wh_proj_sigs(:,2),wh_proj_sigs(:,3));

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
            for m=midxs
                cst = medoid_cost(cl_assns,wh_proj_sigs(:,m),wh_proj_sigs,midxs);
                if cst < min_cost
                    min_cost    = cst;
                    min_m       = m;
                end

            end
            
            size(cl_meds(:,j))
            size(min_m)
            cl_meds(:,j) = wh_proj_sigs(:,min_m);
        end

    end

end

function cst = medoid_cost(cl_assns,medoid,wh_proj_sigs,midxs)

    cst = 0;
    for i=midxs

        cst = cst + norm(wh_proj_sigs(:,i)-medoid);
        
    end
        
end