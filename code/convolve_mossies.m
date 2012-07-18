function PsVm = convolve_mossies(GC_model,mean_mf)

time=size(mean_mf,2);
padtime=200;

mf_pad=[mean_mf(:,1:padtime) mean_mf(:,1:padtime) mean_mf];

Ps=convolve_matrix_by_tau(GC_model.dt,GC_model.tau_s,mf_pad);

%how about this?
% slowrectifying=convolve_matrix_by_tau(GC_model.dt,30,mf_pad);
% %Ps=(Ps+slowrectifying)/2;
% rat=GC_model.tau_s_ratio;
% Ps=rat*Ps+(1-rat)*slowrectifying;


PsVm=convolve_matrix_by_tau(GC_model.dt,GC_model.tau_m,Ps);

PsVm=PsVm(:,1+padtime*2:time+padtime*2);