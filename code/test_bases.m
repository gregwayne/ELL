% [celltypes, bfns]=generate_bases('raw');
% [celltypes, bfns]=generate_bases('2thr',.5, .8); %first is PCA, 2nd is others
[celltypes, bfns]=generate_bases('thr',.8);
ncells=size(bfns,1);

smk=500;
smk2=1000;
rates=bfns;
indother=cellfun(@isempty,strfind(celltypes,'PCA'));
indPCA=setdiff(1:ncells,find(indother));
indPCAstrict=find(strcmp(celltypes,'PCA'));
% rates(indPCA,400:900)=rates(indPCA,400:900)*5;

smk=10000;
smk2=5000;
rates(indPCA,:)=smoothts(rates(indPCA,:),'g',smk2*5,smk2);
rates(indother,:)=smoothts(rates(indother,:),'g',smk*5,smk);
% rates=rates+.05;


% rates(indPCA,:)=[];ncells=size(rates,1);

vars=zeros(ncells,1);
vars=mean(rates,2);
% vars(strcmp(celltypes,'PCA'))=2;
% vars(~strcmp(celltypes,'PCA'))=50;

NCs=zeros(ncells);

pCI=.1; %pretend some cells have strong common input!
indCI=randperm(ncells);%indPCA(randperm(length(indPCA)));
indCI=indCI(1:length(indCI)*pCI);
NCs(indCI,:)=1;
NCs=NCs.*NCs';


C=zeros(ncells);
% C=C+NCs;
% C=C+diag(vars);

% rtemp=rates;
% rtemp(indother,:)=rtemp(indPCA,:)/20;
% rtemp=rates-rtemp;
% C=rtemp*rtemp'+diag(vars);

% C=rates*rates'+NCs+diag(vars);

k=makekernel(rates,C);
% k=k-ones(4500,1)*mean(k,1);



%plot stuff!
tran=(1:4500)*5e-5;
figure(1);imagesc(tran,1:124,rates);title('basis functions!')
figure(2);imagesc(tran,tran,k);title('kernel!')
check_LR(k);