MFearly=~cellfun(@isempty,strfind(mftypes,'early'));
MFmed=~cellfun(@isempty,strfind(mftypes,'med'));
MFlate=~cellfun(@isempty,strfind(mftypes,'late'));
MFpause=~cellfun(@isempty,strfind(mftypes,'pause'));

Wuse=Wsparse;

early=find(sum(bsxfun(@times,Wuse ,MFearly)'~=0));
med=find(sum(bsxfun(@times,Wuse,MFmed)'~=0));
late=find(sum(bsxfun(@times,Wuse,MFlate)'~=0));
pause=find(sum(bsxfun(@times,Wuse,MFpause)'~=0));

Enn=length(setdiff(early,union(union(med,late),pause)));
EM=length(setdiff(intersect(early,med),union(late,pause)));
EL=length(setdiff(intersect(early,late),union(med,pause)));
EP=length(setdiff(intersect(early,pause),union(med,late)));

Mnn=length(setdiff(med,union(union(early,late),pause)));
ML=length(setdiff(intersect(med,late),union(early,pause)));
MP=length(setdiff(intersect(med,pause),union(early,late)));

Lnn=length(setdiff(late,union(union(med,early),pause)));
LP=length(setdiff(intersect(late,pause),union(early,med)));

Pnn=length(setdiff(pause,union(union(early,med),late)));

EML=length(setdiff(intersect(intersect(early,med),late),pause));
EMP=length(setdiff(intersect(intersect(early,med),pause),late));
ELP=length(setdiff(intersect(intersect(early,late),pause),med));
MLP=length(setdiff(intersect(intersect(late,med),pause),early));


nE      = 99;
nM      = 8;
nL      = 2;
nP      = 9;
nEM     = 16;
nEL     = 11;
nEP     = 15;
nML     = 0;
nLP     = 0;
nPM     = 0;
nEMP    = 3;
nELM    = 0;
nELP    = 0;
nLMP    = 0;
nN      = 12; %gonna disregard these

all_terms_meas      = [nE,nM,nL,nP,nEM,nEL,nEP,nML,nPM,nLP,nELM,nEMP,nELP,nLMP];
all_terms_fit       = [Enn,Mnn,Lnn,Pnn,EM,EL,EP,ML,MP,LP,EML,EMP,ELP,MLP];

all_terms_meas=all_terms_meas/sum(all_terms_meas);
all_terms_fit=all_terms_fit/sum(all_terms_fit);

bardata=[all_terms_meas;all_terms_fit];
clf;bar(bardata')
set(gca,'XTickLabel',{'Early', 'Medium', 'Late', 'Pause', 'E+M', 'E+L', 'E+P', 'M+L', 'P+M', 'L+P', 'E+L+M', 'E+M+P', 'E+L+P', 'L+M+P'})
legend('Observed','In model fits')


