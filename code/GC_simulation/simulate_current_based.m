function [realtrace,modeltrace,tran] = simulate_current_based(GC_model,mean_mf,real_cells)

dt          = GC_model.dt;
min_t       = GC_model.min_t;
max_t       = GC_model.max_t;
times       = min_t:dt:max_t;
nsteps      = length(times)-1;
V_gr        = GC_model.E_L*ones(1,nsteps);
Ps          = zeros(3,1);
inputs=zeros(3,nsteps);
trialnum=zeros(3,1);
                  

vmean=zeros(size(V_gr));
for rep=1:GC_model.nreps
    Ps = zeros(3,1);    %initialize synaptic current to 0 at start of each trial
    for i=1:3 %pick another set of trials
        if(GC_model.MF_input(i)>0)
%             trialnum(i)=ceil(rand*size(rspstore{GC_model.MF_input(i)},1));
%             inputs(i,:)=rspstore{GC_model.MF_input(i)}(trialnum(i),1:nsteps);
        inputs(i,:)=mean_mf(GC_model.MF_input(i),:);                                %changed to just use the mean MF input
        end
    end
    for step=1:nsteps-1

        if GC_model.spiking_on
            if V_gr(step) > GC_model.V_thresh
                V_gr(step) = GC_model.V_reset; 
            end
        end

        % increment presynaptic conductance when there's a spike
        Ps = (1 - (dt/GC_model.tau_s))*Ps + inputs(:,step)*dt/GC_model.tau_s;
        
        %V_gr(step+1)  = V_gr(step) + (dt/GC_model.tau_m)*(GC_model.E_L - V_gr(step) - GC_model.Ws*Ps*V_gr(step));
        V_gr(step+1)  = V_gr(step) + (dt/GC_model.tau_m)*(GC_model.E_L - V_gr(step) + GC_model.Ws*Ps);

        if GC_model.spiking_on
            if V_gr(step+1) > GC_model.V_thresh
                V_gr(step) = 0;
            end
        end

    end
    vmean=vmean+V_gr/GC_model.nreps;
end

realtrace=real_cells(GC_model.GC_to_model,1:nsteps);

%use convolution instead of simulation:
% for i=1:3
%     if(GC_model.MF_input(i)>0)
%         inputs(i,:)=mean_mf(GC_model.MF_input(i),:);
%     end
% end
% Ps=convolve_matrix_by_tau(GC_model.dt,GC_model.tau_s,inputs);
% vmean=GC_model.Ws*convolve_matrix_by_tau(GC_model.dt,GC_model.tau_m,Ps);

modeltrace=vmean(1:nsteps);
tran=(min_t+dt:dt:max_t)*10^-3;








