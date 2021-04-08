function data = readtable_array(fname)
%my own function to read table from csv and convert it to array, with
%missing values set as nan

data=table2array(readtable(fname,'ReadVariableNames',true));
if iscellstr(data)
    data = cellfun(@(str) strrep(str,',','.'), data,'UniformOutput',false);
    data = cellfun(@str2num, data,'UniformOutput',false);
    WoR_delay_empty = cell2mat(cellfun(@isempty, data,'UniformOutput',false));
    data(WoR_delay_empty) = {NaN};
    data = cell2mat(data);
end

end