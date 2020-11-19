addpath \\lexport\iss01.charpier\analyses\lgi1\Git-Paul\Persyst_to_Brainvision\PersystImportExport1.0\PersystImportExport1.0

datafile = 'C:\Users\paul.baudin\Pictures\Diapos LGI1\Schwa\test_export.lay';

%fieldtrip cannot read those data
hdr = ft_read_header(datafile);

[hdr, data] = layread(datafile);