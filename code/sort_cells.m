function sorted=sort_cells(data,basis)
%function sorted=sort_cells(data,basis)
%returns cells in data sorted by their projection onto basis

ncells=size(data,1);
ordering=zeros(ncells,1);
for i=1:ncells
    ordering(i)=mean((data(i,:)-mean(data(i,:))).*(basis-mean(basis)));%/(var(data(i,:))*var(basis)));
end
sorted=[ordering data];
sorted=sortrows(sorted,1);
sorted=sorted(:,2:end);