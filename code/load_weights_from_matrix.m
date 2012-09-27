function GC_model=load_weights_from_matrix(GC_model,Wmat,gcnum)

GC_model.Ws=nonzeros(Wmat(gcnum,:));

if(size(GC_model.Ws,1)~=1) %get dimensions right
    GC_model.Ws=GC_model.Ws';
end

GC_model.Ws=padarray(GC_model.Ws,[0 3-length(GC_model.Ws)],'post');
GC_model.MF_input=find(Wmat(gcnum,:));
GC_model.MF_input=padarray(GC_model.MF_input,[0 3-length(GC_model.MF_input)],'post');
GC_model.GC_to_model=gcnum;