function GC_resetweights(source,eventdata,uictrls)
    for i=1:3
        set(uictrls.synapse_sliders{i},'Value',0.5);
    end