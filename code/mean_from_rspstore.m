function meanV=mean_from_rspstore(rspstore)

for i=1:length(rspstore)
    meanV(i,:)=mean(rspstore{i});
end