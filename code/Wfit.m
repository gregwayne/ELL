function GC_model=Wfit(GC_model,useweight,countcap,rspstore,real_cells,v_fig,h1,h2)
%this was a first attempt at fitting weights in the model GCs-- but
%optWeights does better
if(length(useweight)>=1)

    step=4;
    stepsize=1e-2;
    GC_model.Ws(useweight)=.1;
    GC_model=Wfit(GC_model,useweight(2:end),3,rspstore,real_cells,v_fig,h1,h2);
    [realtrace,modeltrace,~] = redoplot(GC_model,rspstore,real_cells);
    GC_model.Ws(useweight)=.1*sqrt(var(realtrace)/var(modeltrace));
    fitold=1000;
    fitnew=fitold-1;
    
    count=1;
    while (step>1&&count<countcap)
        while fitnew<fitold
            step=step*2;
            fitold=fitnew;
            GC_model.Ws(useweight)=GC_model.Ws(useweight)+stepsize*step;
            
            GC_model=Wfit(GC_model,useweight(2:end),3,rspstore,real_cells,v_fig,h1,h2);
            
            [realtrace,modeltrace,tran] = redoplot(GC_model,rspstore,real_cells);
            realtrace=realtrace-mean(realtrace);
            modeltrace=modeltrace-mean(modeltrace);
            fitnew=mean(sqrt((realtrace-modeltrace).^2));
            
            set(0,'CurrentFigure',v_fig);
            set(h1,'YData',realtrace-mean(realtrace(1:200)));
            set(h2,'YData',modeltrace-mean(modeltrace(1:200)));
            title('Adjusting weights...');
            drawnow;
            
        end
        GC_model.Ws(useweight)=GC_model.Ws(useweight)-stepsize*step;
        step=step/2;
        GC_model.Ws(useweight)=GC_model.Ws(useweight)-stepsize*step;
        step=step/4;
        fitold=1000;
        fitnew=fitold-1;
        count=count+1;
    end
    step=step*4;
    GC_model.Ws(useweight)=GC_model.Ws(useweight)+stepsize*step;
    GC_model=Wfit(GC_model,useweight(2:end),3,rspstore,real_cells,v_fig,h1,h2);
end