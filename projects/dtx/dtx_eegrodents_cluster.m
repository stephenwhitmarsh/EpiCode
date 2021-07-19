function dtx_eegrodents_cluster(ipatient, config_script)

%% Set parameters
if ispc
    addpath \\lexport\iss01.charpier\analyses\lgi1\Git-Paul\fieldtrip;
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\projects\dtx'));
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\external'));
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\shared'));
elseif isunix
    addpath /network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/fieldtrip/
    addpath(genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/projects/dtx'));
    addpath(genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/external'));
    addpath(genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/shared'));
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

config = eval(config_script);
ipart = 1;

%% align LFP, remove artefacted trials, correct baseline

MuseStruct = readMuseMarkers(config{ipatient},false);
MuseStruct = dtx_remove_wrong_seizure(config{ipatient}, MuseStruct, false);

config{ipatient}.LFP.name   = {'SlowWave_EMG_begin','Crise_End'};
LFP_others                   = readLFP(config{ipatient},MuseStruct,false);

config{ipatient}.align.name = {'SlowWave'};
config{ipatient}.LFP.name   = {'SlowWave', 'Seizure'};
MuseStruct                  = alignMuseMarkersPeaks(config{ipatient},MuseStruct,false);
LFP_SlowWave                = readLFP(config{ipatient},MuseStruct,false);

LFP{1}.SlowWave             = LFP_SlowWave{1}.SlowWave;
LFP{1}.Seizure              = LFP_SlowWave{1}.Seizure;
LFP{1}.SlowWave_EMG_begin   = LFP_others{1}.SlowWave_EMG_begin;
LFP{1}.Crise_End            = LFP_others{1}.Crise_End;
clear LFP_SlowWave* LFP_others

%flip data if required
for markername = string(config{ipatient}.name)
    if isempty(LFP{ipart}.(markername))
        continue
    end
    if istrue(config{ipatient}.LFP.flip)
        for itrial = 1:size(LFP{ipart}.(markername).trial,2)
            LFP{ipart}.(markername).trial{itrial} = LFP{ipart}.(markername).trial{itrial} .* -1;
        end
    end
end

save(fullfile(config{ipatient}.datasavedir,[config{ipatient}.prefix,'LFP.mat']), 'LFP', '-v7.3');

%% compute slow wave morpho
LFP = removeArtefactedTrials(config{ipatient}, LFP);
datasavedir = config{ipatient}.datasavedir;
markername = "SlowWave";

%lp filter
cfgtemp                   = [];
cfgtemp.lpfilter          = 'yes';
cfgtemp.lpfreq            = 30;
cfgtemp.lpfilttype        = 'fir';
LFP{ipart}.(markername)   = ft_preprocessing(cfgtemp,LFP{ipart}.(markername));

fig = figure;hold on;
for itrial = 1:size(LFP{ipart}.(markername).trial,2)
    cfgtemp = [];
    cfgtemp.trials = itrial;
    LFP_onetrial = ft_selectdata(cfgtemp, LFP{ipart}.(markername));
    
    config{ipatient}.morpho.channame = config{ipatient}.align.channel.(markername);
    try
        [hw,amp] = plot_morpho(config{ipatient}, LFP_onetrial);
    catch
        warning('cannot find slowwave in trial %d', itrial);
        hw = nan;
        amp = nan;
    end
    
    cfgtemp         = [];
    cfgtemp.latency = config{ipatient}.LFP.baselinewindow.(markername);
    cfgtemp.channel = config{ipatient}.align.channel.(markername);
    data_sel        = ft_selectdata(cfgtemp, LFP_onetrial);
    std_bl          = std(data_sel.trial{:}, 'omitnan');
    
    morpho.halfwidth(itrial) = hw;
    morpho.amplitude(itrial) = amp;
    morpho.amplitude_norm(itrial) = amp/std_bl;
    startdir                 = MuseStruct{1}{LFP_onetrial.trialinfo.idir}.starttime;
    morpho.time(itrial)      = startdir + seconds((LFP_onetrial.trialinfo.begsample - LFP_onetrial.trialinfo.offset) / LFP_onetrial.fsample);
end
%remove text to make the figure readable
delete(findall(gcf,'type','text'));

%print to file
fname = fullfile(config{ipatient}.imagesavedir,'..','morpho_eachrat',sprintf('%smorpho_%s',config{ipatient}.prefix,markername));
dtx_savefigure(fig,fname,'pdf','png','close');

%save computed morpho
fprintf('save morpho values to %s\n',fullfile(datasavedir,sprintf('%sslowwave_morpho.mat',config{ipatient}.prefix)));
save(fullfile(datasavedir,sprintf('%sslowwave_morpho.mat',config{ipatient}.prefix)), 'morpho', '-v7.3');




