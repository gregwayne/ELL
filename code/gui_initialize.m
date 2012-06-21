clear all;
close all;

%load your data (run this once)
[mftypes,rspstore]=export_mossy_rasters('../final_mossyfibers/');
[gctypes,real_cells]=generate_bases('raw');