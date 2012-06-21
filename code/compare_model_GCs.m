%load your data (run this once)
[mftypes,rspstore]=export_mossy_rasters('../final_mossyfibers/');
[gctypes,real_cells]=generate_bases('raw');
%%

tau_m       = 8.7; % 8.7 (ms)
V_thresh    = -43; % (mV);
rmgs        = 1; % unitless
E_L         = -63; % (mV)
V_reset     = E_L; % (mV);
tau_s       = 0.25; % (ms) -- look up
N_m         = 3;
Ws          = [1.2705 0 0];
Ws=Ws/sqrt(tau_s);
incells = [21 0 0]; %0 for no input
cell_to_model=22;

dt          = 5e-2; % (ms)
min_t       = -0.025e3;
max_t       = .2e3;
times       = min_t:dt:max_t;
nsteps      = length(times)-1;
V_gr        = E_L*ones(1,nsteps);
Ps          = zeros(N_m,1);    
spiking_on  = 0;
inputs=zeros(N_m,nsteps);
trialnum=zeros(N_m,1);
                  

vmean=zeros(size(V_gr));
nreps=10;
for rep=1:nreps
    for i=1:N_m %pick another set of trials
        if(incells(i)>0)
            trialnum(i)=ceil(rand*size(rspstore{incells(i)},1));
            inputs(i,:)=rspstore{incells(i)}(trialnum(i),1:nsteps);
        end
    end
    for step=1:nsteps-1

        if spiking_on
            if V_gr(step) > V_thresh
                V_gr(step) = V_reset; 
            end
        end

        % increment presynaptic conductance when there's a spike
        Ps              = Ps + inputs(:,step); 
        Ps              = Ps - (dt/tau_s)*Ps;

        V_gr(step+1)  = V_gr(step) + (dt/tau_m)*(E_L - V_gr(step) - Ws*Ps*V_gr(step));

        if spiking_on
            if V_gr(step+1) > V_thresh
                V_gr(step) = 0;
            end
        end

    end
    vmean=vmean+V_gr/nreps;
end

clf;hold on;
plot((min_t+dt:dt:max_t)*10^-3,real_cells(cell_to_model,1:nsteps)-mean(real_cells(cell_to_model,1:200)))
plot(.5e-3+(min_t+dt:dt:max_t)*10^-3,vmean-mean(vmean(1:200)),'g');
legend('real cell','fake cell','location','northeast')
xlim([-.025 .05])
fixfig