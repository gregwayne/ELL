function GC_model=fitterLS(GC_model,mean_mf,real_cells)

H=convolve_mossies(GC_model,mean_mf);
inputs=H(nonzeros(GC_model.MF_input),:);
v=real_cells(GC_model.GC_to_model,:);

%add some offset inputs
inputs=[inputs;ones(1,length(inputs))];
inputs=[inputs;-1*ones(1,length(inputs))];

% w=v*pinv(inputs);
w=lsqnonneg(inputs',v');

count=1;
for i=1:3
    if(GC_model.MF_input(i)~=0)
        GC_model.Ws(i)=w(count);
        count=count+1;
    else
        GC_model.Ws(i)=0;
    end
end