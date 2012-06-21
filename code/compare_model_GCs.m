function redoplot(GC_model,rspstore,real_cells)

dt          = 5e-2; % (ms)
min_t       = -0.025e3;
max_t       = .2e3;
times       = min_t:dt:max_t;
nsteps      = length(times)-1;
V_gr        = GC_model.E_L*ones(1,nsteps);
Ps          = zeros(3,1);    
inputs=zeros(3,nsteps);
trialnum=zeros(3,1);
                  

vmean=zeros(size(V_gr));
nreps=10;
for rep=1:nreps
    for i=1:3 %pick another set of trials
        if(incells(i)>0)
            trialnum(i)=ceil(rand*size(rspstore{incells(i)},1));
            inputs(i,:)=rspstore{incells(i)}(trialnum(i),1:nsteps);
        end
    end
    for step=1:nsteps-1

        if spiking_on
            if V_gr(step) > GC_model.V_thresh
                V_gr(step) = GC_model.V_reset; 
            end
        end

        % increment presynaptic conductance when there's a spike
        Ps              = Ps + inputs(:,step); 
        Ps              = Ps - (dt/GC_model.tau_s)*Ps;

        V_gr(step+1)  = V_gr(step) + (dt/GC_model.tau_m)*(GC_model.E_L - V_gr(step) - GC_model.Ws*Ps*V_gr(step));

        if spiking_on
            if V_gr(step+1) > GC_model.V_thresh
                V_gr(step) = 0;
            end
        end

    end
    vmean=vmean+V_gr/nreps;
end

cla;hold on;
plot((min_t+dt:dt:max_t)*10^-3,real_cells(GC_model.GC_to_model,1:nsteps)-mean(real_cells(GC_model.GC_to_model,1:200)))
plot(.5e-3+(min_t+dt:dt:max_t)*10^-3,vmean-mean(vmean(1:200)),'g');
legend('real cell','fake cell','location','northeast')
xlim([-.025 .15])
fixfig