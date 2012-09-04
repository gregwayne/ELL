function input = draw_MF_input(mossyraster)

trnum = randperm(size(mossyraster,1));
input = mossyraster(trnum(1),:);