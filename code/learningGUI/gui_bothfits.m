function gui_bothfits(source,eventdata,uictrls,mean_mf,real_cells,mftypes)

gui_fitweights_lasso(source,eventdata,uictrls,mean_mf,real_cells,mftypes);
drawnow;
gui_fitweights_LS(source,eventdata,uictrls,mean_mf,real_cells,mftypes);