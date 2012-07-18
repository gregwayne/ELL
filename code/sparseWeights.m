function [GC_model,dropCount]=sparseWeights(Wstore,usecell,mean_mf,real_cells,bonus)

GC_model_initialize();
useweights=find(Wstore(usecell,:));
inputs=convolve_mossies(GC_model,mean_mf(useweights,:));
v=real_cells(usecell,:);

%deal our weights
ncells=length(useweights);
maxwts=6; %do this better blah
wtable=[];
for i=ncells:-1:1
    dat=nchoosek(1:ncells,i);
    if(size(dat,2)~=i)
        dat=dat';
    end
    wtable=[wtable; padarray(dat,[0 maxwts-i],'post')];
end

%run LS on all combinations
mse=zeros(size(wtable,1),1);
inputs_shift=[];
for i=1:size(wtable,1)
    inputs_shift=[inputs(nonzeros(wtable(i,:)),:);ones(1,length(inputs))];
    inputs_shift=[inputs_shift;-1*ones(1,length(inputs))];
    [lsfits,mse(i)]=lsqnonneg(inputs_shift',v');
    w(i,1:maxwts)=padarray(lsfits(1:end-2)',[0 maxwts-length(lsfits)+2],'post');
end

%and now compare
mse_adjusted=mse - mse(1)*((1+bonus)*(ncells-sum(wtable~=0,2)));
[~,ind]=min(mse_adjusted);

GC_model.GC_to_model=usecell;
GC_model.MF_input=padarray(useweights(nonzeros(wtable(ind,:))),[0 maxwts-sum(wtable(ind,:)~=0)],'post');
GC_model.Ws=w(ind,:);

dropCount=sum(wtable(1,:)~=0)-sum(wtable(ind,:)~=0);