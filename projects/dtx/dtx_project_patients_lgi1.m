
%% Analyse of LGI1 patients seizures
%Paul Baudin

%% Set parameters
addpath C:\Users\paul.baudin\Documents\MATLAB\fieldtrip;
addpath (genpath('C:\Users\paul.baudin\Documents\MATLAB\EpiCode\projects\dtx'));
addpath (genpath('C:\Users\paul.baudin\Documents\MATLAB\EpiCode\external'));
addpath (genpath('C:\Users\paul.baudin\Documents\MATLAB\EpiCode\shared'));
addpath C:\Users\paul.baudin\Documents\MATLAB\DTX;
addpath C:\Users\paul.baudin\Documents\MATLAB\MatlabImportExport_v6.0.0;
ft_defaults
feature('DefaultCharacterSet', 'CP1252') % To fix bug for weird character problems in reading neuralynx

[config] = dtx_setparams_patients_lgi1([]);


for ipatient = 1:length(config)
    
    
    %% Get LFP data
    
    [MuseStruct]               = readMuseMarkers(config{ipatient}, false);
    [MuseStruct]               = alignMuseMarkers(config{ipatient},MuseStruct, false);
    [dat_LFP]                  = readLFP(config{ipatient}, MuseStruct, false, false);
    %dat_LFP{ipart}{imarker}
    
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
    
    
    [dat_LFP, MuseStruct] = dtx_removeartefactsLFP(dat_LFP, MuseStruct, config{ipatient}, false);
    %remove dat_LFP et dat_EMG : ce sera data{imarker}
    
    
    %% Plot
    
    
    
    for ipart = 1:length(dat_LFP)
        
        fprintf('*********************************************\n');
        fprintf('*********************************************\n');
        fprintf('**** Plot data, for %s, part %d ****\n',config{ipatient}.prefix(1:end-1), ipart);
        fprintf('*********************************************\n');
        fprintf('*********************************************\n\n');
        
        
        dataEEG = []; %for figure with right and left eeg data, without emg data
        iEEG = 0;
            
        %1 figure per marker type
        for imarker = 1:length(config{ipatient}.LFP.name)
            
            
            if ~isempty(dat_LFP{ipart}{imarker})
                
                if any(strfind(config{ipatient}.LFP.name{imarker}, 'EMG')) %is EMG
                    
                    dtx_plot_emg_method(config{ipatient},dat_LFP,ipart,imarker,true);
                    
                    dtx_plot_timecourse_eeg_emg(config{ipatient}, dat_LFP, ipart, imarker, 'emg', true);
                    
                    [EEGwithEMG_onechan_avg{ipatient}{imarker}, EMG_onechan__avg{ipatient}{imarker}] = ...
                        dtx_plot_comparison_eeg_emg(config{ipatient}, dat_LFP, ipart, imarker, true);
                    
                else %is EEG
                    
                    dtx_plot_timecourse_eeg_emg(config{ipatient}, dat_LFP, ipart, imarker, 'eeg', true);
                    
                    dtx_plot_overdraw_allchannels(config{ipatient},dat_LFP, ipart, imarker, true);
                    
                    [EEG_avg_allchan{ipatient}{imarker}, EEG_align_avg{ipatient}{imarker}]  =...
                        dtx_plot_avg_allchannels(config{ipatient},dat_LFP,ipart,imarker,true);
                    
                    dtx_plot_SlowWaveTopographyTimecouse(config{ipatient},dat_LFP, ipart, imarker, true);
                    
                    %[TFR_eeg{ipatient}] = dtx_plot_TFR(config{ipatient}, dat_LFP, ipart, imarker, true);
                    
                    iEEG = iEEG + 1;
                    dataEEG{iEEG} = dat_LFP{ipart}{imarker};
                    dataEEG{iEEG}.LFP.name = config{ipatient}.LFP.name;
                    
                end
                
              
                
            end
        end %imarker
        
        %1 figure common for all markers (SW_R and SW_L)
        
        dtx_plot_SlowWaveTopography(config{ipatient},dataEEG,ipart, true);
        
        if config{ipatient}.continous == true
            Seizure_Infos{ipatient} =...
                dtx_plot_count_seizure(config{ipatient}, MuseStruct, ipart, true);
        else
            Seizure_Infos{ipatient} = [];
        end
                
    end %ipart
    clear dat_LFP dataEEG MuseStruct
end %ipatient


%% save data all patients

%save. Same datasavedir for all patients.
%data are the data of the last ipart : so the 'merged' data
save(fullfile(config{1}.datasavedir,'All_patients_EEG_avg_allchan.mat'),'EEG_avg_allchan');
save(fullfile(config{1}.datasavedir,'All_patients_EEG_align_avg.mat'),'EEG_align_avg');
save(fullfile(config{1}.datasavedir,'All_patients_EMG_onechan__avg.mat'),'EMG_onechan__avg');
save(fullfile(config{1}.datasavedir,'All_patients_EEGwithEMG_onechan_avg'),'EEGwithEMG_onechan_avg');
save(fullfile(config{1}.datasavedir,'All_patients_Seizure_Infos.mat'),'Seizure_Infos');


%% Figure all patients

[config] = dtx_setparams_patients_lgi1([]);

% Same datasavedir for all patients
load(fullfile(config{1}.datasavedir,'All_patients_EEG_avg_allchan.mat'),'EEG_avg_allchan');
load(fullfile(config{1}.datasavedir,'All_patients_EEG_align_avg.mat'),'EEG_align_avg');
load(fullfile(config{1}.datasavedir,'All_patients_EMG_onechan__avg.mat'),'EMG_onechan__avg');
load(fullfile(config{1}.datasavedir,'All_patients_EEGwithEMG_onechan_avg'),'EEGwithEMG_onechan_avg');
load(fullfile(config{1}.datasavedir,'All_patients_Seizure_Infos.mat'),'Seizure_Infos');

%attendre d'avoir une première version des fichiers sauvegardés
%appenddata et voir si c'est OK, si les avg se sont bien concaténés en
%différents trials
%Puis appliquer les fonctions précédentes
%peut être normaliser toutes les traces

%Figure comparaison EEG EMG (pas forcément R et L)
%Figure timecourse all chans avg : R et L
%Figure timecourse chan align (pas forcément R et L)
%Figure timecourse chan EMG (pas forcément R et L)
%Figure topo : R et L + évolution topo (+ film ?)


