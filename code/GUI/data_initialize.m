%this code loads all of the mossy fiber and GC data from the mat files Nate
%gave us (run it once to initialize your workspace before using other
%functions)
%
%   mftypes is a cell array of strings giving names of the mossy fibers
%
%   gctypes is a cell array of strings giving names of the granule cells
%
%   rspstore is a cell array of matrices which stores the raw MF data.
%   Each cell corresponds to a mossy fiber (rspstore{1} is mftypes{1},
%   etc), and contains a matrix of spiking data with N rows (corresponding
%   to trials) and T columns (corresponding to time, using Nate's
%   typical window of -0.025s : 5e-5s : 0.200s relative to the EOCD.)
%
%   real_cells is an N by T matrix of average granule cell responses, in
%   which N is the number of granule cells and T is time (over the same
%   window as above).
%
%   numGCs and numMFs are just what it says on the box

[mftypes,rspstore]=loader_mossies('../mossyfibers_final/','alltrials');
[gctypes,real_cells]=generate_bases('raw','../gcs_mat');

mean_mf=zeros(length(rspstore),size(rspstore{1},2));
for i=1:length(rspstore)
    mean_mf(i,:)=mean(rspstore{i});
end


numGCs=length(gctypes);
numMFs=length(mftypes);

%scale mossies by the dt of nate's data-- so that a trial with one spike
%has integral 1.
data_dt = 0.05;
mean_mf=mean_mf/data_dt;
for i=1:length(rspstore)
    rspstore{i}=rspstore{i}/data_dt;
end


% we increased the time resolution of simulation to make convolution more
% accurate, so better start adjusting all our loaded data at bootup so we
% don't make mistakes later on.


rspstore={}; %too big
[real_cells,mean_mf,rspstore] = interp_data_timestep(real_cells,mean_mf,rspstore,10);

%scratch that, interpolating rspstore takes up too much memory. figure out
%something more efficient!