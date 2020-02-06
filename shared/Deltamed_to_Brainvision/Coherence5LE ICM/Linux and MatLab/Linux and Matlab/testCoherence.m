clear all;
close all;

addpath(genpath('\\lexport\iss01.epimicro\rodents\raw\Scripts_Deltamed_to_vhdr'));
datafile = '\\lexport\iss01.epimicro\rodents\raw\Scripts_Deltamed_to_vhdr\Coherence5LE ICM\Data files\8ch-4096Hz.EEG'
z = Coherence5LE()
z = OpenFile(z, datafile);

% eeg_folder = 'C:\Users\thomas.bancel\Documents\2018_matlab_thomas_internship\data_deltamed\EEG2\';
% eeg_filename = '180413a-b_0004.eeg';
% z = OpenFile(z,strcat(eeg_folder, eeg_filename));

markers = GetMarkersNumber(z);
disp([num2str(markers) ' Markers found']);
[z,res] = PutMarker(z,100,0,'test',1);
markers = GetMarkersNumber(z);
disp(['Putmarker returned: ' num2str(res)]);
disp([num2str(markers) ' Markers found']);
FileInfo = GetFileInfo(z);
disp('Getting Markers');
markers = GetMarkers(z,0,FileInfo.frequency*FileInfo.duration)
disp('Getting Impedances');
imps = GetImpedances(z,FileInfo.frequency*FileInfo.duration)

disp('Plotting data ...');
data = GetData(z);
plot(data);

% only works if file is a bloc-file
%disp('Trying NextFile');
%z = NextFile(z,1)