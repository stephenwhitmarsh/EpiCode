% 
addpath /Users/stephen.whitmarsh/WhitmarshEpilepsy/
addpath /Users/stephen.whitmarsh/fieldtrip/ 
ft_defaults

patient_directory       = '/Users/stephen.whitmarsh';
patient_directory       = '/Volumes/epimicro/Donnees-analyses/Stephen/pat_02230_0674/eeg';
directory_searchstring  = '02230_2015-02-25_*';
directory_searchstring  = '02230_2015-02-*';
data_searchstring       = '*m1*.ncs';

MuseStruct = readMuseMarkers(patient_directory,directory_searchstring,data_searchstring);

writeMuseMarkers(MuseStruct);

writeSpykingCircusDeadFile(MuseStruct);