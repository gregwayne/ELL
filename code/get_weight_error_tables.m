function [mse, w, wtable] = get_weight_error_tables(Wstore,usecell,mean_mf,real_cells)

GC_model_initialize();

useweights=find(Wstore(usecell,:));
inputs=convolve_mossies(GC_model,mean_mf(useweights,:));
v=real_cells(usecell,:);
v=v(~isnan(v));

%deal our weights
ncells=length(useweights);
maxwts=3; %cap the number of nonzero weights you can fit
wtable=[];
for i=min(ncells,maxwts):-1:1
    dat=nchoosek(1:ncells,i);
    if(size(dat,2)~=i)
        dat=dat';
    end
    wtable=[wtable; padarray(dat,[0 maxwts-i],'post')];
end

%run LS on all combinations
mse = zeros(size(wtable,1),1);
inputs_shift = [];
w=zeros(size(wtable,1),maxwts);
for i=1:size(wtable,1)
    inputs_shift    = [inputs(nonzeros(wtable(i,:)),:);ones(1,length(inputs))];
    inputs_shift    = [inputs_shift;-1*ones(1,length(inputs))];
    inputs_shift    = inputs_shift(:,1:length(v));
    [lsfits,~] = lsqnonneg(inputs_shift',v');
    mse(i) = var(v - lsfits'*inputs_shift)/var(v);
    w(i,1:maxwts)   = padarray(lsfits(1:end-2)',[0 maxwts-length(lsfits)+2],'post');
end

