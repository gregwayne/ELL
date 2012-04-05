% [celltypes, bfns]=generate_bases('raw');
[celltypes, bfns]=generate_bases('thr',.8);
ncells=size(bfns,1);

smk=10000;
rates=smoothts(bfns,'g',smk*5,smk);


vars=zeros(ncells,1);
% vars=mean(rates,2);
% vars(strcmp(celltypes,'PCA'))=1;
% vars(~strcmp(celltypes,'PCA'))=50;

NCs=zeros(ncells);

C=zeros(ncells);

rtemp=rates;
rtemp(cellfun(@isempty,strfind(celltypes,'PCA')),:)=0;
% rtemp=rates-rtemp;
C=rtemp*rtemp'/1e9;

% C=rates*rates'+NCs+diag(vars);

k=makekernel(rates,C);
% k=k-ones(4500,1)*mean(k,1);



%plot stuff!
figure(1);imagesc(rates);title('basis functions!')
% figure(2);imagesc(k);title('kernel!')
check_LR(k);