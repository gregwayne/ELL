function file_list = get_file_list(generic_filenames)

num_files = 0;
file_list = {};
for i=1:length(generic_filenames)
    fn = char(generic_filenames(i));
    [folder, ~, ~] = fileparts(fn);

    filenames = dir(fn);
    for j=1:length(filenames)
        num_files = num_files + 1;
        file_list(num_files) = {strcat(folder, '/', filenames(j).name)};
    end
end

end