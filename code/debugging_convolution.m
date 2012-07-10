
GC_model=get_model_from_gui(uictrls)

mfconv=convolve_matrix_by_tau(GC_model.dt,GC_model.tau_s,mean_mf);
H=convolve_matrix_by_tau(GC_model.dt,GC_model.tau_m,mfconv);
H=H(:,1:4500);

dat=real_cells(GC_model.GC_to_model,:);
H=bsxfun(@minus,H,mean(H,2));
[w,stats]=lasso(H',dat,'DFmax',3);
[w2,stats2]=lasso(H',dat,'DFmax',20);
w=w(:,1);
w2=w2(:,1);
%%
clc
disp('weights fit with DFmax=3');
[nonzeros(w) find(w)]
disp(' ');
disp('weights fit with DFmax=20');
[nonzeros(w2) find(w2)]

%%
figure(5);clf;hold on;
tran=GC_model.min_t+GC_model.dt:GC_model.dt:GC_model.max_t;
plot(tran/1e3,dat-mean(dat(1:100)));

GC_model_new=GC_model;
nonzerocells=find(w);
newWs=w(nonzerocells);
for i=1:3
    GC_model_new.Ws(i)=newWs(i);
    GC_model_new.MF_input(i)=nonzerocells(i);
end

[~,modeltrace,~,~] = current_based(GC_model,mean_mf,real_cells);
plot(tran/1e3,modeltrace-mean(modeltrace(1:100)),'r');

[~,modeltrace2,~,~] = current_based(GC_model_new,mean_mf,real_cells);
plot(tran/1e3,modeltrace2-mean(modeltrace2(1:100)),'c--');

fitcell=H'*w;
plot(tran/1e3,fitcell-mean(fitcell(1:100)),'g');

fitcell2=H'*w2;
plot(tran/1e3,fitcell2-mean(fitcell2(1:100)),'m');

% ylim([-2 14])
xlim([GC_model.min_t GC_model.max_t]/1e3)
legend('real cell','fit in gui, euler','fit outside of gui, euler','fit outside of gui, convolution','fit outside, DFmax=20')
fixfig