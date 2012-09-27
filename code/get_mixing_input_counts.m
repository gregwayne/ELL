function counts = get_mixing_input_counts(W,mftypes,numcats)

MF_indices = get_struct_of_celltypes(mftypes);

early   = sum(W(:,MF_indices.early),2)~=0;
med     = sum(W(:,MF_indices.med),2)~=0;
late    = sum(W(:,MF_indices.late),2)~=0;
pause   = sum(W(:,MF_indices.pause),2)~=0;


nE = sum(early.*(1-med).*(1-late).*(1-pause));
nM = sum(med.*(1-early).*(1-late).*(1-pause));
nL = sum(late.*(1-early).*(1-med).*(1-pause));
nP = sum(pause.*(1-early).*(1-med).*(1-late));

nEM = sum(early.*med.*(1-late).*(1-pause));
nEL = sum(early.*late.*(1-med).*(1-pause));
nEP = sum(early.*pause.*(1-med).*(1-late));
nML = sum(med.*late.*(1-early).*(1-pause));
nPM = sum(pause.*med.*(1-early).*(1-late));
nLP = sum(late.*pause.*(1-early).*(1-med));

nELM = sum(early.*late.*med.*(1-pause));
nEMP = sum(early.*med.*pause.*(1-late));
nELP = sum(early.*late.*pause.*(1-med));
nLMP = sum(late.*med.*pause.*(1-early));

if(numcats==4)
    counts = [nE,nM,nL,nP,nEM,nEL,nEP,nML,nPM,nLP,nELM,nEMP,nELP,nLMP];
else
    counts = [nE,nM+nL+nML,nP,nEM+nEL+nELM,nEP,nPM+nLP+nLMP,nEMP+nELP];
end