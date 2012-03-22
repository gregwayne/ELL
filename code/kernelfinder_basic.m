pth='../DATA2';
rcell=loader_means(pth);
tmax=3000;
rmat=cell2mat(cellfun(@(x) [x(1:min(length(x),tmax))' x(end)*ones(1,max(tmax-length(x),0))]',rcell,'UniformOutput',false))';

ncells=size(rmat,1);

rmat=rmat(randperm(size(rmat,1)),:);
rbar=rmat(1:ncells,:);
C=rbar*rbar';
Cinv=C^-1;
K=rbar'*Cinv*rbar;
D=zeros(tmax*2-1,1);
for i=1:tmax*2-1
    D(i)=sum(diag(K,i-tmax));
end


% let's take a look!
figure(1);
subplot(2,1,1);
imagesc(K)
ylabel('t1')
xlabel('t2')
subplot(2,1,2);
plot(-tmax+1:tmax-1,D);
xlabel('Offset')
ylabel('Sum(diag(K))')