function Wfit = balanced(GC_model)

    tau_s = GC_model.tau_s;
    tau_m = GC_model.tau_m;
    
    ts = 0:GC_model.dt:10;
    
    V_bar = -GC_model.E_L;
    E_L   = GC_model.E_L;
    
    %wo    = 1;

    %Wfit  = wo*V_bar;
    lms     = log(tau_m/tau_s); 
    
    e1      = exp(-lms*tau_m/(tau_m-tau_s));
    e2      = exp(-lms*tau_s/(tau_m-tau_s));
    te      = (1/(tau_s-tau_m))*(e1-e2);
    Wfit    = 1/te;
    
    wo      = Wfit/V_bar;
    
    Vs = bexps(ts,tau_m,tau_s,E_L,V_bar,wo);
%     plot(ts,Vs);
    
    t_max = log(tau_m/tau_s)*(tau_s*tau_m)/(tau_m-tau_s);
%     
%     disp(sprintf('t_max is %f',t_max));
%     disp(sprintf('Wfit is %f',Wfit));

end

function y = bexps(t,tau_m,tau_s,E_L,V_bar,w)

    y   = E_L + w*V_bar*(1/(tau_s-tau_m))*(exp(-t/tau_s) - exp(-t/tau_m));
    
end