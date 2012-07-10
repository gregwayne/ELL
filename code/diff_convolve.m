function y = diff_convolve(v,tau,dt)

ts=0:dt:100;
kernel=exp(-ts/tau)/tau;
y=zeros(100/dt,1);

for i=1:length(v)
    y(i+1)=(1-(dt/tau))*y(i) + (dt/tau)*v(i);
end