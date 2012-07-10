function GC_model=gridfitter(GC_model,mean_mf,real_cells,wmin,wmax,wabsvals,gridres,nreps)

for rep=1:nreps
    mincost=inf;
    nmf=length(nonzeros(GC_model.MF_input));
    usept=zeros(nmf,1);
    gridspan = wmax(1)-wmin(1);

    grid=makeGrid(nmf,gridres,wmin,wmax);
    for i=1:length(grid)
        count=1;
        for j=1:3 %deal the weights
            if(GC_model.MF_input(j)>0)
                GC_model.Ws(j)=grid(i,count);
                count=count+1;
            end
        end

        %compute error at current gridpoint using MSE
        costnew=compute_model_error(GC_model,mean_mf,real_cells,'MSE');

        if(costnew<mincost)
            usept=grid(i,:);
            mincost=costnew;
        end
    end

    count=1;
    for j=1:3
        if(GC_model.MF_input(j)>0)
            GC_model.Ws(j)=usept(count);
            count=count+1;
        end
    end

    if(nreps>1)
        wmin=usept-gridspan/gridres;
        wmax=usept+gridspan/gridres;
        wmin(wmin<min(wabsvals))=min(wabsvals);
        wmax(wmax>max(wabsvals))=max(wabsvals);
    end
end