function grid=makeGrid(ndim,res,min,max)
%function grid=makeGrid(ndim,res,min,max)
% ndim = number of input cells
% res = grid resolution
% min, max = ndim x 1 vectors giving grid range for each cell
grid=[];
if ndim>1
    subgrid=makeGrid(ndim-1,res,min(1:ndim-1),max(1:ndim-1));
else
    grid=linspace(min(ndim),max(ndim),res)';return;
end
for x=min(ndim):(max(ndim)-min(ndim))/(res-1):max(ndim)
    for i=1:size(subgrid,1)
        grid=[grid; x subgrid(i,:)];
    end
end