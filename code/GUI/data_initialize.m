[mftypes,rspstore]=export_mossy_rasters('../mossyfibers_final/');
[gctypes,real_cells]=generate_bases('raw');
mean_mf=mean_from_rspstore(rspstore);