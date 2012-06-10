granule_cell_filenames = {};
for i=1:length(gc_names)
   
    name = gc_names(i);
    granule_cell_filenames = {granule_cell_filenames{:},'../DATA2/PCA/*.mat'};
   
end

gc_files    = get_file_list(granule_cell_filenames);
gcs         = {};
E_L         = -70; % (mV)

%%
for i=1:length(gc_files)
    
   gcs{i} = importdata(char(gc_files(i)));
   start  = gcs{i}.start;
   inter  = gcs{i}.interval;
   leng   = gcs{i}.length;
   times  = start:inter:(leng*inter + start - inter);
 
   plot(times,gcs{i}.values);
   xlim([0 times(end)]);
   ylim([E_L -30]);
   pause(1);
   
end

%% Pull out closest match in Nate's database
%% make the parameter files be separate