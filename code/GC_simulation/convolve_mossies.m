function PsVm = convolve_mossies(GC_model,mean_mf)

time=size(mean_mf,2);
padtime=2000;

mf_pad=[mean_mf(:,1:padtime) mean_mf(:,1:padtime) mean_mf];

Ps=convolve_matrix_by_tau(GC_model.dt,GC_model.tau_s,mf_pad);
PsVm=convolve_matrix_by_tau(GC_model.dt,GC_model.tau_m,Ps);

PsVm=PsVm(:,1+padtime*2:time+padtime*2);