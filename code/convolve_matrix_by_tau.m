function convolved_matrix = convolve_matrix_by_tau(dt,tau,input_matrix)
%takes dt,tau,input_matix, where each row of input_matrix is a spike
%train/time series. Convolves rows of input_matrix with an exponential
%kernel with time constant tau; does not attenuate the result, so returned
%matrix is longer than original.

tran=0:dt:100;
kernel=(1/tau)*exp(-tran/tau);

%add padding at the front?
convolved_matrix=zeros(size(input_matrix,1),size(input_matrix,2)+length(tran)-1);
for i=1:size(input_matrix,1)
    convolved_matrix(i,:)=conv(input_matrix(i,:),kernel)*dt;
end