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

isPatient = false;
isEEGvideo = false;
isProbe = false;

isPatient = false;

if slurm_task_id<=3
    ispeak = true;
else
    ispeak = false;
end

if slurm_task_id == 1 || slurm_task_id == 4
    config = dtx_setparams_patients_lgi1([]);
    if ispeak
        for i=1:length(config)
            config{i}.name{1} = 'SlowWave_R_peak';
            config{i}.name{2} = 'SlowWave_L_peak';
            config{i}.LFP.name = config{i}.name;
        end
    end
    isPatient = true;
elseif slurm_task_id == 2 || slurm_task_id == 5
    config = dtx_setparams_eegvideo([]);
    if ispeak
        for i=1:length(config)
            config{i}.name{1} = 'SlowWave_peak';
            config{i}.LFP.name = config{i}.name;
        end
    end
    isEEGvideo = true;
elseif slurm_task_id == 3 || slurm_task_id == 6
    [config] = dtx_setparams_probe_lfp([]);
    isProbe = true;
    if ispeak
        for i=1:length(config)
            config{i}.name{1} = 'SlowWave_EEG_peak';
            config{i}.LFP.name = config{i}.name;
        end
    end
end

for ipatient = 1:length(config)
    %% TEMPORAIRE PAUL
    config{ipatient}.imagesavedir = fullfile(config{ipatient}.imagesavedir,config{ipatient}.name{1});
    config{ipatient}.datasavedir = fullfile(config{ipatient}.datasavedir, config{ipatient}.name{1});
    %
    %% Get LFP data
    
    [MuseStruct]                    = readMuseMarkers(config{ipatient}, true);
    MuseStruct_origin               = MuseStruct; %for counting seizures, it does not have to with the seizures removed by miss alignment
    [MuseStruct]                    = alignMuseMarkers_peak(config{ipatient},MuseStruct, true);
    [MuseStruct]                    = alignMuseMarkers_begin(config{ipatient},MuseStruct,true);
    [MuseStruct]                    = alignMuseMarkers_EMG(config{ipatient},MuseStruct,true);
    [dat_LFP]                       = readLFP(config{ipatient}, MuseStruct, true, true); %dat_LFP{ipart}{imarker}
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
    [dat_LFP, MuseStruct] = removetrialLFP_marker(config{ipatient}, dat_LFP, MuseStruct,'all', 'all','remove','BAD__START__','BAD__END__',0,0,0,0, true);
    
    % correct baseline
    for ipart = 1:size(dat_LFP,2)
        for imarker = 1:size(dat_LFP{ipart},2)
            if ~isempty(dat_LFP{ipart}{imarker})
                cfgtemp = [];
                cfgtemp.demean = 'yes';
                cfgtemp.baselinewindow = config{ipatient}.align.toibaseline{1};
                dat_LFP{ipart}{imarker} = ft_preprocessing(cfgtemp, dat_LFP{ipart}{imarker});
            end
        end
    end
    %% Plot patient by patient
    
    
    
    for ipart = 1:length(dat_LFP)
        
        
        dataEEG = []; %for figure with right and left eeg data, without emg data
        iEEG = 0;
        
        %1 figure per marker type
        for imarker = 1:length(config{ipatient}.LFP.name)
            
            if ~isempty(dat_LFP{ipart}{imarker}) %some imarker do not have data in some patients
                
                %SlowWavealign : EMG data aligned from SlowWave. Remove all
                %trials without EMG
                if any(strfind(config{ipatient}.LFP.name{imarker}, 'SlowWavealign'))
                    [dat_LFP, MuseStruct] = removetrialLFP_marker(config{ipatient}, dat_LFP, MuseStruct,ipart,imarker,'keep',config{ipatient}.LFP.emgmarker{imarker},config{ipatient}.LFP.emgmarker{imarker},0,0,-5,5, true);
                end
                
                
                
                if ~isempty(dat_LFP{ipart}{imarker}.trial)
                    
                    electrodeToPlot{1} = config{ipatient}.LFP.electrodetoplot{imarker};
                    if isPatient
                        electrodeToPlot{2} = config{ipatient}.LFP.motorcortex{imarker};
                    end
                    
                    for i=1:length(electrodeToPlot)
                        dtx_plot_timecourse_eeg_emg(config{ipatient}, dat_LFP, ipart, imarker, 'eeg',electrodeToPlot{i},config{ipatient}.epoch.toi{imarker}, true);
                        dtx_plot_timecourse_eeg_emg(config{ipatient}, dat_LFP, ipart, imarker, 'eeg',electrodeToPlot{i},[-5, 5], true);
                    end
                    
                    if isfield (config{ipatient}.LFP, 'TFR')
                        if config{ipatient}.LFP.TFR.doTFR == true
                            dtx_plot_TFR(config{ipatient}, dat_LFP, ipart, imarker, false, true);
                        end
                    end
                    
                    if isPatient
                        dtx_plot_SlowWaveTopographyTimecouse(config{ipatient},dat_LFP, ipart, imarker, true);
                    end
                    
                    dtx_plot_avg_allchannels(config{ipatient},dat_LFP,ipart,imarker,config{ipatient}.epoch.toi{imarker},true);
                    dtx_plot_avg_allchannels(config{ipatient},dat_LFP,ipart,imarker,[-2 2],true);
                    dtx_plot_overdraw_allchannels(config{ipatient},dat_LFP, ipart, imarker,config{ipatient}.epoch.toi{imarker}, true);
                    dtx_plot_overdraw_allchannels(config{ipatient},dat_LFP, ipart, imarker,[-2 2], true);
                    %keep eeg for plot of both SW_R and SW_L in the same figure
                    iEEG = iEEG + 1;
                    dataEEG{iEEG} = dat_LFP{ipart}{imarker};
                    dataEEG{iEEG}.LFP.name = config{ipatient}.LFP.name{imarker};
                    
                    
                    %plot only imarkers associated with EMG
                    if any(strfind(config{ipatient}.LFP.name{imarker}, 'EMG'))
                        
                        %choose electrodes
                        emgToPlot = config{ipatient}.LFP.emg{imarker};
                        if isPatient
                            eegToPlot = config{ipatient}.LFP.motorcortex{imarker};
                        else
                            eegToPlot = config{ipatient}.LFP.electrodetoplot{imarker};
                        end
                        
                        %do plots
                        dtx_plot_emg_method(config{ipatient},dat_LFP,ipart,imarker,emgToPlot,true);
                        dtx_plot_timecourse_eeg_emg(config{ipatient}, dat_LFP, ipart, imarker, 'emg',emgToPlot,config{ipatient}.EMG.toi, true);
                        dtx_plot_comparison_eeg_emg(config{ipatient}, dat_LFP, ipart, imarker, true, eegToPlot, emgToPlot);
                    end
                    
                    
                    
                    %get data for pool all patients
                    [data_avg_allchans{imarker}{ipatient}, data_avg_chantoplot{imarker}{ipatient}, data_TFR{imarker}{ipatient}, data_EMG{imarker}{ipatient}] =...
                        dtx_get_patient_data(config{ipatient}, dat_LFP, ipart, imarker);
                    
                    
                else
                    data_avg_allchans{imarker}{ipatient}        = [];
                    data_avg_chantoplot{imarker}{ipatient}      = [];
                    data_TFR{imarker}{ipatient}                 = [];
                    data_EMG{imarker}{ipatient}                 = [];
                end
            else
                data_avg_allchans{imarker}{ipatient}        = [];
                data_avg_chantoplot{imarker}{ipatient}      = [];
                data_TFR{imarker}{ipatient}                 = [];
                data_EMG{imarker}{ipatient}                 = [];
            end %if dat_LFP{ipart}{imarker} is empty
            
        end %imarker
        
        %Figure common for all markers
        
        if isPatient
            if length(dataEEG) >= 2
                dtx_plot_SlowWaveTopography(config{ipatient},dataEEG(1:2),ipart, true); %only the 2 firsts 'dataeeg' for patients
            else
                dtx_plot_SlowWaveTopography(config{ipatient},dataEEG,ipart, true);
            end
        end
       
        if config{ipatient}.continous == true %Based on manually-put marker, leave markers with artefacts, and with alignment non detected
            Seizure_Infos{ipatient} = dtx_plot_count_seizure(config{ipatient}, MuseStruct_origin, ipart, false);
        else
            Seizure_Infos{ipatient} = [];
        end
        
        
        datalist = [];
        channames = [];
        datalist = dat_LFP{ipart};
        channames = config{ipatient}.LFP.electrodetoplot(1:length(datalist));
        dtx_plot_comparison_severaleeg(config{ipatient},datalist,ipart,channames, [-2 2], true, true)
        
        
    end %ipart
    clear dat_LFP dataEEG MuseStruct
end %ipatient


%% save data all patients

%save. Same datasavedir for all patients.
%data are the data of the last ipart : so the 'merged' data
save(fullfile(config{1}.datasavedir,'All_patients_data_avg_allchansn.mat'),'data_avg_allchans');
save(fullfile(config{1}.datasavedir,'All_patients_data_avg_chantoplot.mat'),'data_avg_chantoplot');
save(fullfile(config{1}.datasavedir,'All_patients_data_TFR.mat'),'data_TFR');
save(fullfile(config{1}.datasavedir,'All_patients_data_EMG.mat'),'data_EMG');
save(fullfile(config{1}.datasavedir,'All_patients_Seizure_Infos.mat'),'Seizure_Infos');



%% Figures all patients
% A REVOIR

% keep same data structure as for patient by patient analysis
ipart       = 1;
ipatient    = 1;

%set config
config{ipatient}.prefix             = 'AllPatients-';
config{ipatient}.imagesavedir       = fullfile(config{ipatient}.imagesavedir,'..','..','AllPatients',config{ipatient}.name{1});
config{ipatient}.merge              = false;

%rename all EMG to plot as 'EMG'
if isfield(config{ipatient}.LFP, 'emg')
    for i = 1:length(config{ipatient}.LFP.emg)
        if ~strcmp(config{ipatient}.LFP.emg{i},'no')
            config{ipatient}.LFP.emg{i} = 'EMG';
        end
    end
end


% Append data : one patient is one trial. Ignore empty 'imarker'
for imarker = 1:length(config{ipatient}.LFP.name)
    
    isdata = [];
    isdata = ~cellfun('isempty',data_avg_allchans{imarker}); %remove empty arrays for append
    if any(isdata)
        data_avg_allchans_allpatients{ipart}{imarker}                = ft_appenddata([],data_avg_allchans{imarker}{isdata});
    else
        data_avg_allchans_allpatients{ipart}{imarker}                = [];
    end
    
    isdata = [];
    isdata = ~cellfun('isempty',data_avg_chantoplot{imarker}); %remove empty arrays for append
    if any(isdata)
        data_avg_chantoplot_allpatients{ipart}{imarker}                  = ft_appenddata([],data_avg_chantoplot{imarker}{isdata});
    else
        data_avg_chantoplot_allpatients{ipart}{imarker}                = [];
    end
    
    isdata = [];
    isdata = ~cellfun('isempty',data_EMG{imarker}); %remove empty arrays for append
    if any(isdata)
        data_EMG_allpatients{ipart}{imarker}        = ft_appenddata([],data_EMG{imarker}{isdata});
    else
        data_EMG_allpatients{ipart}{imarker}                = [];
    end
    
end

% for do_normalize = [false true]
%normalize by dividing by value of max SlowWave at t=0 (alignment).
%normalize emg by dividing by max(EMG) beetween 0 and 2
%     if do_normalize
%         config{ipatient}.prefix = 'AllPatients-normalized-';
%         config{ipatient}.imagesavedir = [config{ipatient}.imagesavedir,'-normalized'];
%
%         for imarker = 1:length(config{ipatient}.LFP.name)
%
%             for itrial = 1:length(config) %one trial is one patient
%                 %set norm values
%                 t = [];
%                 norm_eeg = 0;
%                 norm_emg = 0;
%                 if ~isempty(data_avg_chantoplot_allpatients{ipart}{imarker})
%                     if itrial <= length(data_avg_chantoplot_allpatients{ipart}{imarker}.trial)
%                         t_debut = data_avg_chantoplot_allpatients{ipart}{imarker}.time{itrial}>-0.5;
%                         t_fin   = data_avg_chantoplot_allpatients{ipart}{imarker}.time{itrial}<0.5;
%                         t_norm  = logical(t_debut .* t_fin);
%                         norm_eeg = max(data_avg_chantoplot_allpatients{ipart}{imarker}.trial{itrial}(t_norm));
%                     end
%                 end
%                 if ~isempty(data_EMG_allpatients{ipart}{imarker})
%                     if itrial <= length(data_EMG_allpatients{ipart}{imarker}.trial)
%                         motorcortex_indx = find(strcmp(data_EMG_allpatients{ipart}{imarker}.label, config{ipatient}.LFP.motorcortex{imarker}));
%                         t = data_EMG_allpatients{ipart}{imarker}.time{itrial};
%                         norm_emg = max(data_EMG_allpatients{ipart}{imarker}.trial{itrial}(end,t>=-2 & t<=2)); %end = EMG, defined by dtx_get_patient_data
%                     end
%                 end
%
%                 %apply norm values to data
%                 if ~isempty(data_avg_allchans_allpatients{ipart}{imarker})
%                     if itrial <= length(data_avg_allchans_allpatients{ipart}{imarker}.trial)
%                         data_avg_allchans_allpatients{ipart}{imarker}.trial{itrial} = data_avg_allchans_allpatients{ipart}{imarker}.trial{itrial} / norm_eeg;
%                     end
%                 end
%
%                 if ~isempty(data_EMG_allpatients{ipart}{imarker})
%                     if itrial <= length(data_EMG_allpatients{ipart}{imarker}.trial)
%                         data_EMG_allpatients{ipart}{imarker}.trial{itrial}(motorcortex_indx,:) = data_EMG_allpatients{ipart}{imarker}.trial{itrial}(motorcortex_indx,:) / norm_eeg;
%                         data_EMG_allpatients{ipart}{imarker}.trial{itrial}(end,:) = data_EMG_allpatients{ipart}{imarker}.trial{itrial}(end,:) / norm_emg;
%                     end
%                 end
%
%                 if ~isempty(data_avg_chantoplot_allpatients{ipart}{imarker})
%                     if itrial <= length(data_avg_chantoplot_allpatients{ipart}{imarker}.trial)
%                         data_avg_chantoplot_allpatients{ipart}{imarker}.trial{itrial} = data_avg_chantoplot_allpatients{ipart}{imarker}.trial{itrial} / norm_eeg;
%                     end
%                 end
%
%             end
%         end %imarker
%
%     end %do_normalize


%plot data : exactly the same protocol as patient by patient
dataEEG = []; %for figure with right and left eeg data, without emg data
iEEG = 0;

for imarker = 1:length(config{ipatient}.LFP.name)
    
    %EMG analysis
    if ~isempty(data_EMG_allpatients{ipart}{imarker})
        
        emgToPlot = 'EMG';
        if isPatient
            eegToPlot = config{ipatient}.LFP.motorcortex{imarker};
        else
            eegToPlot = config{ipatient}.LFP.electrodetoplot{imarker};
        end
        dtx_plot_emg_method(config{ipatient},data_EMG_allpatients,ipart,imarker,emgToPlot,true);
        dtx_plot_timecourse_eeg_emg(config{ipatient}, data_EMG_allpatients, ipart, imarker, 'emg',emgToPlot,config{ipatient}.EMG.toi, true);
        dtx_plot_comparison_eeg_emg(config{ipatient}, data_EMG_allpatients, ipart, imarker, true, eegToPlot, emgToPlot);
        dtx_plot_comparison_eeg_emg_trialbytrial(config{ipatient}, data_EMG_allpatients, ipart, imarker, true, eegToPlot, emgToPlot)
        
    end
    
    %EEG analysis
    
    if ~isempty(data_avg_allchans_allpatients{ipart}{imarker})
        
        electrodeToPlot = config{ipatient}.LFP.name{imarker};
        dtx_plot_timecourse_eeg_emg(config{ipatient}, data_avg_chantoplot_allpatients, ipart, imarker, 'eeg',electrodeToPlot,config{ipatient}.epoch.toi{imarker}, true);
        dtx_plot_timecourse_eeg_emg(config{ipatient}, data_avg_chantoplot_allpatients, ipart, imarker, 'eeg',electrodeToPlot,[-5 5], true);
        
        if isPatient
            electrodeToPlot = config{ipatient}.LFP.motorcortex{imarker};
            dtx_plot_timecourse_eeg_emg(config{ipatient}, data_avg_allchans_allpatients, ipart, imarker, 'eeg',electrodeToPlot,config{ipatient}.epoch.toi{imarker}, true);
            dtx_plot_timecourse_eeg_emg(config{ipatient}, data_avg_allchans_allpatients, ipart, imarker, 'eeg',electrodeToPlot,[-5 5], true);
        end
        
        
        if isfield (config{ipatient}.LFP, 'TFR')
            if config{ipatient}.LFP.TFR.doTFR == true
                dtx_plot_TFR(config{ipatient}, data_TFR, ipart, imarker, true, true);
            end
        end
        
        if isPatient
            dtx_plot_SlowWaveTopographyTimecouse(config{ipatient},data_avg_allchans_allpatients, ipart, imarker, true);
        end
        
        dtx_plot_overdraw_allchannels(config{ipatient},data_avg_allchans_allpatients, ipart, imarker,config{ipatient}.epoch.toi{imarker}, true);
        dtx_plot_overdraw_allchannels(config{ipatient},data_avg_allchans_allpatients, ipart, imarker,[-2 2], true);
        dtx_plot_avg_allchannels(config{ipatient},data_avg_allchans_allpatients,ipart,imarker,config{ipatient}.epoch.toi{imarker},true);
        dtx_plot_avg_allchannels(config{ipatient},data_avg_allchans_allpatients,ipart,imarker,[-2 2],false);
        
        %keep eeg for plot of both SW_R and SW_L in the same figure
        iEEG = iEEG + 1;
        dataEEG{iEEG} = data_avg_allchans_allpatients{ipart}{imarker};
        dataEEG{iEEG}.LFP.name = config{ipatient}.LFP.name{imarker};
        
    end
    
end %imarker

if isPatient
    %1 figure common for all markers (SW_R and SW_L)
    if length(dataEEG) >= 2
        dtx_plot_SlowWaveTopography(config{ipatient},dataEEG(1:2),ipart, true);
    else
        dtx_plot_SlowWaveTopography(config{ipatient},dataEEG,ipart, true);
    end
end

imarker = 1;
%notempty = ~cellfun(@isempty,data_avg_chantoplot_allpatients{ipart}(1:2));
%data_avg_allchanstoplot_allpatients{ipart}{imarker} = ft_appenddata([], data_avg_chantoplot_allpatients{ipart}{notempty});
datalist = [];
channames = [];
datalist = data_avg_chantoplot_allpatients{ipart};
channames = config{ipatient}.LFP.name(1:length(datalist));
dtx_plot_comparison_severaleeg(config{ipatient},datalist,ipart,channames, [-2 2], true, true);


end


% end %do_normalize



