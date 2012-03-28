pth='../DATA2';
rcell=loader_means(pth);
tmax=4000;

rmat=cell2mat(cellfun(@(x) [x(1:min(length(x),tmax))' x(end)*ones(1,max(tmax-length(x),0))],rcell(:,2),'UniformOutput',false));
rmat=rmat-mean(rmat,2)*ones(1,tmax);

ncells=size(rmat,1);
C=rmat*rmat'/tmax^2;

%let's make up some noise correlations!
nc_same=0;
nc_diff=0;
nc_self=0;
NC=zeros(ncells);
for i=1:ncells
    NC(i,i)=nc_self;
    for j=i+1:ncells
        if(strcmp(rcell(i,1),rcell(j,1)))
            NC(i,j)=nc_same; NC(j,i)=nc_same;
        else
            NC(i,j)=nc_diff; NC(j,i)=nc_diff;
        end 
    end
end

C=C+NC;
Cinv=C^-1;
K=rmat'*Cinv*rmat;
D=zeros(tmax*2-1,1);
for i=1:tmax*2-1
    D(i)=sum(diag(K,i-tmax));
end

% let's take a look!
figure(1);
subplot(2,1,1);
imagesc(K)
subplot(2,1,2);
plot(-tmax+1:tmax-1,D);