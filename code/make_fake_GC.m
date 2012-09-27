function GC_model = make_fake_GC(GCparams,MF_indices)
GC_model_initialize;

clear modeltype;
emptycount=3;
while(emptycount>2)
    emptycount=0;
    for i=1:3
        inclass = rand(1);
        if(inclass<GCparams.theta.e)
            modeltype{i} = 'early';
        elseif(inclass<(GCparams.theta.e+GCparams.theta.m))
            modeltype{i} = 'med';
        elseif(inclass<(GCparams.theta.e+GCparams.theta.m+GCparams.theta.l))
            modeltype{i} = 'late';
        elseif(inclass<(GCparams.theta.e+GCparams.theta.m+GCparams.theta.l+GCparams.theta.p))
            modeltype{i} = 'pause';
        else
            modeltype{i} = 'none';
            emptycount=emptycount+1;
        end
    end
end

for i=1:3
    if(~strcmp(modeltype{i},'none'))
        GC_model.MF_input(i) = MF_indices.(modeltype{i})(ceil(rand(1)*length(MF_indices.(modeltype{i}))));
        GC_model.Ws(i) = gamrnd(GCparams.W_gamma_k.(modeltype{i}),GCparams.W_gamma_th.(modeltype{i}));
    else
        GC_model.MF_input(i) = 0;
        GC_model.Ws(i) = 0;
    end
end

GC_model.V_thresh = GC_model.E_L + max(GCparams.thresh_mean + randn(1)*sqrt(GCparams.thresh_var),1e-6); %don't let the threshold go below the resting potential!