
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

[config, mergeindex] = dtx_setparams_patients_lgi1([]);

ipart = 1; %script not adapted yet to data with several "ipart"

for merge_eeg = true %[false true] %one analysis patient per patient, one other merging
    

    if merge_eeg == false
        patients_list = 1:length(config);
    elseif merge_eeg == true
        patients_list = 1:length(mergeindex);
    end
    
    
    for ipatient = patients_list%1:nb_of_patients %go patient/eeg by patient/eg

        
        %% Get LFP data
        
        if merge_eeg == false
            idata = ipatient;
        elseif merge_eeg == true
            idata = mergeindex{ipatient};
        end
        
        i_eeg_to_merge = 0;
        data_temp = [];
        MuseStruct_temp = [];
        
        for ianalyse = idata
            
            force = true;
            if merge_eeg == true
                force = false;
            end
            
            [MuseStruct]    = readMuseMarkers(config{ianalyse}, force);
            [MuseStruct]    = alignMuseMarkers(config{ianalyse},MuseStruct, force);
            [dat_LFP]       = readLFP(config{ianalyse}, MuseStruct, force, force); %dat_LFP{ipart}{imarker}
            
            %inverse the data according to standard clinical visualisation
            for imarker = 1:size(dat_LFP{ipart},2)
                for itrial = 1:size(dat_LFP{ipart}{imarker}.trial,2)
                    dat_LFP{ipart}{imarker}.trial{itrial} = -dat_LFP{ipart}{imarker}.trial{itrial};
                end
            end
            
            if merge_eeg == true
                i_eeg_to_merge = i_eeg_to_merge + 1;
                data_temp{i_eeg_to_merge} = dat_LFP;
                MuseStruct_temp{i_eeg_to_merge} = MuseStruct;
            end
        end
        
        if merge_eeg == true 
            [config{ipatient}, MuseStruct, dat_LFP_EMG, dat_LFP] = dtx_merge_patientslgi1(mergeindex,ipatient,MuseStruct_temp,data_temp);
        else
            dat_LFP_EMG = dat_LFP;
        end
        

        %% Plot
        
        fprintf('*********************************************\n');
        fprintf('*********************************************\n');
        fprintf('**** Plot data, for %s ****\n',config{ipatient}.prefix(1:end-1));
        fprintf('*********************************************\n');
        fprintf('*********************************************\n\n');
        
        
        %1 figure per marker type
        for imarker = 1:length(config{ipatient}.LFP.name)
            
            dtx_plot_timecourse_eeg_emg(config{ipatient}, dat_LFP{ipart},imarker, 'eeg', true);
            dtx_plot_timecourse_eeg_emg(config{ipatient}, dat_LFP_EMG{ipart},imarker, 'emg', true);
            
            [TFR_eeg{ipatient}] = dtx_plot_TFR(config{ipatient}, dat_LFP{ipart}, imarker, true);
            
            dtx_plot_overdraw_allchannels(config{ipatient},dat_LFP{ipart},imarker,true);
            
            EEG_avg_allchan{ipatient}{imarker} = dtx_plot_avg_allchannels(config{ipatient},dat_LFP{ipart},imarker,true);
            
            dtx_plot_SlowWaveTopographyTimecouse(config{ipatient},dat_LFP{ipart}, imarker, true);
            
            if isfield(config{ipatient}.LFP,'emg')
                if ~(strcmp(config{ipatient}.LFP.emg{imarker},'no'))
                    if ~(config{ipatient}.LFP.emg{imarker} == false)
                        
                        [EEGchanalign_avg{ipatient}, EMG_avg{ipatient}] = dtx_plot_comparison_eeg_emg(config{ipatient}, dat_LFP_EMG{ipart}, imarker,true); %find which strategy of EMG quantif
                    
                    end
                end
            end
            
        end
        
        %1 figure common for all markers (SW_R and SW_L)
        
        dtx_plot_SlowWaveTopography(config{ipatient},dat_LFP{ipart},true,true);
        
        if config{ipatient}.continous == true
            Seizure_Infos{ipatient} = dtx_plot_count_seizure(config{ipatient}, MuseStruct, true);
        end
        
    end %ipatient
end %merge_eeg

%% save data all patients

%Remove data of non-merged eeg
config_merged                   = config{1:length(mergeindex));
TFR_eeg                         = TFR_eeg(1:length(mergeindex)); 
EEG_avg_allchan                 = EEG_avg_allchan(1:length(mergeindex)); 
EEGchanalign_avg                = EEGchanalign_avg(1:length(mergeindex)); 
EMG_avg                         = EMG_avg(1:length(mergeindex)); 
Seizure_Infos                   = Seizure_Infos(1:length(mergeindex)); 

%save. Same datasavedir for all patients
save(fullfile(config{1}.datasavedir,'config_merged.mat'),'config_merged');
save(fullfile(config{1}.datasavedir,'All_patients_TFR_eeg.mat'),'TFR_eeg');
save(fullfile(config{1}.datasavedir,'All_patients_EEG_avg_allchan.mat'),'EEG_avg_allchan');
save(fullfile(config{1}.datasavedir,'All_patients_EEGchanalign_avg.mat'),'EEGchanalign_avg');
save(fullfile(config{1}.datasavedir,'All_patients_EMG_avg.mat'),'EMG_avg');
save(fullfile(config{1}.datasavedir,'All_patients_Seizure_Infos.mat'),'Seizure_Infos');


%% Figure all patients

[config, mergeindex] = dtx_setparams_patients_lgi1([]);

% Same datasavedir for all patients
load(fullfile(config{1}.datasavedir,'config_merged.mat'),'config');
load(fullfile(config{1}.datasavedir,'All_patients_TFR_eeg.mat'),'TFR_eeg');
load(fullfile(config{1}.datasavedir,'All_patients_EEG_avg_allchan.mat'),'EEG_avg_allchan'); %récupérer avg de tous les channels, pas un par un
load(fullfile(config{1}.datasavedir,'All_patients_EEGchanalign_avg.mat'),'EEGchanalign_avg');
load(fullfile(config{1}.datasavedir,'All_patients_EMG_avg.mat'),'EMG_avg');
load(fullfile(config{1}.datasavedir,'All_patients_Seizure_Infos.mat'),'Seizure_Infos');


%attendre d'avoir une première version des fichiers sauvegardés
%appenddata et voir si c'est OK, si les avg se sont bien concaténés en
%différents trials
%Puis appliquer les fonctions précédentes