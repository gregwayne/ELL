clear all;
close all;
mossy_fiber_filenames = {   '../mossyfibers/UBC_type1/*.mat'  , 
                            '../mossyfibers/pause/*.mat'      ,
                            '../mossyfibers/pca_early/*.mat'  ,
                            '../mossyfibers/pe_med/*.mat'     };

N_cells = 0;                        
mossy_fibers    = {};
lrange          = 1;  
for i=1:length(mossy_fiber_filenames) 
    new_mossies         = get_file_list(get_file_list({mossy_fiber_filenames{i}}));
    mossy_fibers        = {mossy_fibers{:},new_mossies{:}};
    N_cells             = N_cells + length(new_mossies);
    mossy_types{i}      = [lrange (lrange + length(new_mossies))];
    lrange              = lrange + 1 + length(new_mossies);
end

mossy_recordings        = {};

for i=1:N_cells
   
    m                   = load(mossy_fibers{i});
    mossy_channels{i}   = fieldnames(m);
    mossy_recordings{i} = {};

    for j=1:2
                  
        % put Ch9 or Ch10 in first index, Ch60 in second
        idx = 1 + isempty(findstr('Ch60',char(mossy_channels{i}(j))));  
        mossy_recordings{i}{idx} = getfield(m,mossy_channels{i}{j});           
               
    end
    
end

% Need to deal with recordings of different lengths,
% so divide each cell into individual trials.
% mossy_trials will contain every trial for each cell
mossy_trials = {};
for i=1:N_cells
   
    mossy_trials{i} = {};
    cmd_times       = mossy_recordings{i}{1}.times;
    sp_times        = mossy_recordings{i}{2}.times;
    
    for j=1:(length(cmd_times)-1)
               
        trial_idx = find(sp_times >= cmd_times(j) ...
                        & sp_times <= cmd_times(j+1));
        mossy_trials{i}{j} = sp_times(trial_idx) - cmd_times(j);
               
    end

end

% Theta vector (e,m,p): 
% (0.45,0.06,0.05)

N_gr        = 1e2;
N_m         = 3;
c_table     = zeros(N_gr,N_m);
for i=1:N_gr
    for j=1:N_m
        
        mossy_kind = rand;
        if mossy_kind <= 0.45 % early
           
            crange      = [min(mossy_types{3},mossy_types{4}),...
                            max(mossy_types{3},mossy_types{4})];            
        elseif mossy_kind <= 0.51 % medium
            crange      = mossy_types{1}; 
        elseif mossy_kind <= 0.56 % pause
            crange      = mossy_types{2};
        else % no input
            crange      = [NaN NaN];
        end
                
        %c_table(i,j) = 1+floor(rand*N_cells);
        if isnan(crange(1))
            c_table(i,j)    = NaN;
        else
            c_table(i,j)    = crange(1) + floor(rand*(crange(2)-crange(1)));
        end
    end
end
        
% Granule Parameters
% deterministic for now
tau_m       = 20; % 8.7 (ms)
V_thresh    = -43; % (mV);
rmgs        = 1; % unitless
E_L         = -63; % (mV)
V_reset     = E_L; % (mV);
tau_s       = 1; % (ms) -- look up
Ws          = rmgs*(rand(N_gr,N_m));

N_cycles    = 10;
dt          = 1e-2; % (ms)
max_t       = 200;  % (ms)
times       = 0:dt:max_t;
V_gr        = E_L*ones(N_gr,length(times));
spiking_on  = 0;

% Storage
granule_traces = {};

for cycle=1:N_cycles
    
    disp(sprintf('Cycle # %d', cycle));
    
    % The presynaptic partner for a granule cell
    % is always the same mossy fiber, but the 
    % fiber's information is randomly selected from 
    % all trials for which we have recordings of it
    P       = zeros(N_gr,N_m);
    for i=1:N_gr
       
        for j=1:N_m
           
            idx     = c_table(i,j);
            if isnan(idx)
                partner = NaN;
            else                
                partner = 1+floor(rand*length(mossy_trials{idx}));
            end
            P(i,j)  = partner;
            
        end
        
    end
                 
    inputs  = zeros(N_gr,N_m,length(times));    
    
    % make event times into indices of binary vectors
    % e.g., [1;0;0;0;1;...] means an event on time 1 and 4
    for i=1:N_gr

        for j=1:N_m

            inputs(i,j,:)   = 0;
            if isnan(P(i,j))
                channel         = [];
            else
                channel         = mossy_trials{c_table(i,j)}{P(i,j)};
            end
            channel             = channel*(1000); % convert to (ms) from s
            for event=1:length(channel)
               
                idx             = 1+floor(channel(event)/dt);
                if channel(event) <= 200 % cut off late events
                    inputs(i,j,idx) = 1;
                end
                
            end
            
        end        

    end
        
    % compute granule cell activity       
    for i=1:N_gr
        
        Ps              = zeros(N_m,1);                        
        input_in_steps  = zeros(N_m,length(times));
        for step=1:(length(times)-1)
            
            if spiking_on
                if V_gr(i,step) > V_thresh
                    V_gr(i,step) = V_reset; 
                end
            end
            
            % increment presynaptic conductance when there's a spike
            Ps              = Ps + inputs(i,:,step)'; 
            Ps              = Ps - (dt/tau_s)*Ps;
            
            V_gr(i,step+1)  = V_gr(i,step) ...
                                + (dt/tau_m) ...
                                *(E_L - V_gr(i,step) - Ws(i,:)*Ps*V_gr(i,step));
            
            if spiking_on
                if V_gr(i,step+1) > V_thresh
                    V_gr(i,step) = 0;
                end
            end
        
        end
        
        granule_traces{i,cycle} = V_gr(i,:); 
        
        if 0
            figure(1);
            plot(times,V_gr(i,:),'g');
            xlim([0 max_t]);
            ylim([-70 -30]);
            hold on;
            for j=1:N_m
                scatter(times,V_gr(i,:).*squeeze(inputs(i,j,:))','r','Marker','*');
            end
            hold off;
            pause(1);
        end
        
    end
            
end

sf = '../granule_traces.mat';
save(sf,'granule_traces');