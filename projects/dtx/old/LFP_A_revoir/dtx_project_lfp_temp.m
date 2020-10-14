function dtx_project_lfp(slurm_task_id)
%
% Analysis of LFP data for DTX projects
% Common analysis for patients (patients_lgi1), awake rats (eegvideo) and
% anesthetized rat (probe_lfp)
% Paul Baudin, with the help of Stephen Whitmarsh
%
%
%

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
    addpath /network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/shared
end


ft_defaults

feature('DefaultCharacterSet', 'CP1252') % To fix bug for weird character problems in reading neurlynx


slurm_task_id = 1;

if slurm_task_id == 1
    config = dtx_setparams_patients_lgi1([]);
    isPatient = true;
    
elseif slurm_task_id == 2
    config = dtx_setparams_eegvideo([]);
    config{1}.imagesavedir = fullfile(config{1}.imagesavedir, 'hp_0.1');
    config{1}.datasavedir = fullfile(config{1}.datasavedir, 'hp_0.1');
    config{1}.LFP.hpfilter = 'yes';
    config{1}.LFP.hpfreq = 0.1;
    isEEGvideo = true;
elseif slurm_task_id == 3
    config = dtx_setparams_eegvideo([]);
    config{1}.imagesavedir = fullfile(config{1}.imagesavedir, 'hp_0.5');
    config{1}.datasavedir = fullfile(config{1}.datasavedir, 'hp_0.5');
    config{1}.LFP.hpfilter = 'yes';
    config{1}.LFP.hpfreq = 0.5;
    isEEGvideo = true;
elseif slurm_task_id == 4
    config = dtx_setparams_eegvideo([]);
    config{1}.imagesavedir = fullfile(config{1}.imagesavedir, 'nohp');
    config{1}.datasavedir = fullfile(config{1}.datasavedir, 'nohp');
    config{1}.LFP.hpfilter = 'no';
    isEEGvideo = true;
    
elseif slurm_task_id == 5
    [config] = dtx_setparams_probe_lfp([]);
    isProbe = true;
end


for ipatient = 1:length(config)
    
    
    %% Get LFP data
    
    [MuseStruct]                    = readMuseMarkers(config{ipatient}, false);
    [MuseStruct]                    = alignMuseMarkers(config{ipatient},MuseStruct, false);
    [dat_LFP]                       = readLFP(config{ipatient}, MuseStruct, false, false); %dat_LFP{ipart}{imarker}
    [config{ipatient},dat_LFP]      = dtx_correctDTX2name(config{ipatient},dat_LFP); %correct name of DTX2 (error during acquisition)
    
    
    if config{ipatient}.LFP.flip == true
        %inverse the data according to standard clinical visualisation
        for ipart = 1:size(dat_LFP,2)
            for imarker = 1:size(dat_LFP{ipart},2)
                if ~isempty(dat_LFP{ipart}{imarker})
                    for itrial = 1:size(dat_LFP{ipart}{imarker}.trial,2)
                        dat_LFP{ipart}{imarker}.trial{itrial} = -dat_LFP{ipart}{imarker}.trial{itrial};
                    end
                end
            end
        end
    end
    
    %remove trials with artefacts
    [dat_LFP, MuseStruct] = removetrialLFP_marker(config{ipatient}, dat_LFP, MuseStruct,'all', 'all','remove','BAD__START__','BAD__END__',0,0,0,0, false);
    
    for imarker = 1:2
        halfwidth_SW{ipatient}{imarker} = dtx_plot_SWmorpho(config{ipatient},dat_LFP,length(dat_LFP),imarker,config{ipatient}.LFP.electrodetoplot{imarker}, [-2 2], false, true);
    end
    imarker = 5;
    ipart = length(dat_LFP);
    dtx_plot_comparison_eeg_emg(config{ipatient}, dat_LFP, ipart, imarker, false, 'M1G', config{ipatient}.LFP.emg{imarker});
%     halfwidth_muscle(ipatient) = dtx_plot_SWmorpho(config{ipatient},data_EMG_env,length(dat_LFP),imarker,'EMG', [-2 2], false, true);

end

for imarker = 1:2
    electrodeToPlot = config{ipatient}.LFP.name{imarker};
    dtx_plot_SWmorpho(config{ipatient},data_avg_chantoplot_allpatients,1,imarker,electrodeToPlot, [-2 2], false, true);
end


