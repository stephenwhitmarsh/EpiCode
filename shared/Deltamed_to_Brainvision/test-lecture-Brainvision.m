addpath C:\Users\paul.baudin\Documents\MATLAB\fieldtrip;
addpath (genpath('C:\Users\paul.baudin\Documents\MATLAB\epilepsy'));
addpath C:\Users\paul.baudin\Documents\MATLAB\DTX;
addpath C:\Users\paul.baudin\Documents\MATLAB\MatlabImportExport_v6.0.0;
ft_defaults

feature('DefaultCharacterSet', 'CP1252') % To fix bug for weird character problems in reading neurlynx



data_loaded = ft_read_data('Z:\rodents\raw\DTX-EEGRodents-vhdr\DTXEEG12-100uM\191014r-d_0003.vhdr');
hdr_loaded = ft_read_header('Z:\rodents\raw\DTX-EEGRodents-vhdr\DTXEEG12-100uM\191014r-d_0003.vhdr');

data_loaded = data_loaded*-1;


figure;
hold;
h=10000
for i=1:15
    plot(0:1/hdr_loaded.Fs:20,data_loaded(i,160*4096:180*4096)+h*i);
end

tick = h;
yticks(0 : tick : length(hdr_loaded.label)*h);
set(gca, 'YTickLabel',[hdr_loaded.label(end)', hdr_loaded.label(1:end-1)']); %why Y axis not in the good order
set(gca,'TickDir','out');

