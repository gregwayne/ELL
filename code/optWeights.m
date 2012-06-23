function GC_model=optWeights(uictrls,GC_model,rspstore,real_cells,plotwin,wmin,wmax,gridres,nreps)
mincost=inf;
nmf=length(nonzeros(GC_model.MF_input));
usept=zeros(nmf,1);
gridspan = wmax(1)-wmin(1);

%let the user know what's going on
set(0,'CurrentFigure',plotwin.v_fig);
title(['Adjusting weights...' num2str(nreps)]);
drawnow;

grid=makeGrid(nmf,gridres,wmin,wmax);
for i=1:length(grid)
    count=1;
    for j=1:3 %deal the weights
        if(GC_model.MF_input(j)>0)
            GC_model.Ws(j)=grid(i,count);
            count=count+1;
        end
    end
    
    %compute cost at current point on grid
    [realtrace,modeltrace,~] = redoplot(GC_model,rspstore,real_cells);
    realtrace=realtrace-mean(realtrace);
    modeltrace=modeltrace-mean(modeltrace);
    costnew=mean(sqrt((realtrace-modeltrace).^2)); %using mean squared error
    
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

%comment these lines to suppress visualization of the fitting proces
[realtrace,modeltrace,~] = redoplot(GC_model,rspstore,real_cells);
realtrace=realtrace-mean(realtrace);
modeltrace=modeltrace-mean(modeltrace);
set(plotwin.h1,'YData',realtrace-mean(realtrace(1:200)));
set(plotwin.h2,'YData',modeltrace-mean(modeltrace(1:200)));
drawnow;

if(nreps>1)
    wmin=usept-gridspan/gridres;
    wmax=usept+gridspan/gridres;
    wmin(wmin<0)=get(uictrls.synapse_sliders{1},'Min'); %hacks hacks hacks
    wmax(wmax>10)=get(uictrls.synapse_sliders{1},'Max');
    GC_model=optWeights(uictrls,GC_model,rspstore,real_cells,plotwin,wmin,wmax,gridres,nreps-1);
end