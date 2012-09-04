function cell_indices = get_struct_of_celltypes(celltypes)

cleantypes = regexprep(celltypes,{'\d','\s'},'');
typeinds = unique(cleantypes);

cell_indices = struct;

for i=1:length(typeinds) %logical indexing of cell type membership
    cell_indices.(typeinds{i}) = strcmp(cleantypes,typeinds{i});
end
