function dtx_plot_patients_seizure

% Set parameters
if ispc
    addpath \\lexport\iss01.charpier\analyses\lgi1\Git-Paul\fieldtrip;
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\projects\dtx'));
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\external'));
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\shared'));
elseif isunix
    addpath /network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/fieldtrip/
    addpath(genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/projects/dtx'));
    addpath(genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/external'));
    addpath /network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/shared
end

ft_defaults

%remove fieldtrip's output
ft_warning off
ft_notice off
ft_info off
ft_debug off
global ft_default
ft_default.checkconfig = 'silent';
ft_default.checkpath = 'once';
ft_default.showcallinfo = 'no';
ft_default.trackcallinfo = 'no';
ft_default.tracktimeinfo = 'no';

config      = dtx_patients_lgi1_setparams;

% Example Fig. 1 : 
ipatient    = [11 11];
itrial      = [6 4];
idir        = [6 6];
markername  = ["SlowWave_L", "SlowWave_R"];
ipart       = [1 1];
h           = [50 50];

toi         = [-5 5];
    
%read data
i=1;
MuseStruct = alignMuseMarkersPeaks(config{ipatient(i)}, [], false);
config{ipatient(i)}.LFP.lpfilter   = 'no';
config{ipatient(i)}.LFP.write      = false;
data = readLFP(config{ipatient(i)}, MuseStruct, true);
clear temp

%plot
for i = 1:2
    cfgtemp             = [];
    cfgtemp.lpfilter    = 'yes';
    cfgtemp.lpfreq      = 70;
    cfgtemp.bsfilter    = 'yes';
    cfgtemp.bsfreq      = [49 51];
    data_plot           = ft_preprocessing(cfgtemp, data{ipart(i)}.(markername(i)));
    trialidx = find(data_plot.trialinfo.idir == idir(i) & data_plot.trialinfo.trialnr == itrial(i));
    chanlist = config{ipatient(i)}.LFP.channel; 

    fig = figure; hold on
    color = 'k';
    ichan = 0;
    %eeg
    for channame = string(chanlist)
        if ~any(strcmp(data_plot.label, channame))
            continue
        end
        ichan = ichan+1;
        cfgtemp = [];
        cfgtemp.channel = char(channame);
        data_1chan = ft_selectdata(cfgtemp,data_plot);
        if isempty(data_1chan.label)
            continue
        end
        plot(data_1chan.time{trialidx},data_1chan.trial{trialidx}+(numel(data_plot.label)+1)*h(i)-h(i)*ichan,'Color', color); %first on top
    end
    %emg
    cfgtemp = [];
    cfgtemp.channel = config{ipatient(i)}.EMG.(markername(i));
    EMG = ft_selectdata(cfgtemp,data_plot);
    if ~isempty(EMG.label)
        ichan = ichan+2;
        rescale = 2 * h(i) / (max(EMG.trial{trialidx}) - min(EMG.trial{trialidx}));
        plot(EMG.time{trialidx}, EMG.trial{trialidx}.*rescale + (numel(data_plot.label)+1)*h(i)-h(i)*ichan, 'Color',color); %first on top
    end
    axis tight
    ylim([-h(i) (numel(data_plot.label)+2)*h(i)]);
    xlim(toi);
    xlabel('Time (s)','Interpreter','none', 'Fontsize',15);
    ylabel('Channel name', 'Fontsize',15);
    tick = h(i);
    yticks(0 : tick : (numel(chanlist)+1)*h(i));
    chan_name_all = string([chanlist, {"", 'EMG'}]);
    set(gca, 'FontWeight','bold', 'Fontsize',15, 'FontName', 'Arial', 'YTickLabel',flip(chan_name_all),'TickDir','out');
    title(sprintf('%s : %s, trial %d dir %d',config{ipatient(i)}.prefix(1:end-1),markername(i), itrial(i), idir(i)),'Interpreter','none','Fontsize',18);
    
    fname = fullfile(config{ipatient(i)}.imagesavedir, '..', 'figure_seizure',sprintf('%s%s_dir_%d_trial_%d',config{ipatient(i)}.prefix,markername(i), idir(i), itrial(i)));
    dtx_savefigure(fig, fname, 'png', 'pdf', 'close');
end