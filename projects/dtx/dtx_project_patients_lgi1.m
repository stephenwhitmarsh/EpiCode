
%% Analyse of LGI1 patients seizures
%Paul Baudin
%alignement sur l'électrode dont le pic est le plus clair.
%rejet d'artefact : choix manuel
%EEG : toutes les slowwave sans artefacts
%EMG : seulement les SlowWave avec EMG disponible
%comparaison EEG-EMG : comparaison avec C3 ou C4, le cortex moteur

%à faire : alignement par rapport au début de l'EMG ?
%refaire analyse patient par patient : décommenter ce qui est collé à
%gauche

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

%ajouter la possibilité de calculer sans faire aucun plot

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
                
                %is EMG
%                 if any(strfind(config{ipatient}.LFP.name{imarker}, 'EMG'))
                    
%                     emgToPlot = config{ipatient}.LFP.emg{imarker};               
%                     
%                     dtx_plot_emg_method(config{ipatient},dat_LFP,ipart,imarker,emgToPlot,true);
%                     
%                     dtx_plot_timecourse_eeg_emg(config{ipatient}, dat_LFP, ipart, imarker, 'emg',emgToPlot,config{ipatient}.epoch.toi{imarker}, true);
%                     
%                     dtx_plot_comparison_eeg_emg(config{ipatient}, dat_LFP, ipart, imarker, true, config{ipatient}.LFP.motorcortex{imarker},emgToPlot);
%                     
                    %is EEG
%                 else
                                        
%                     dtx_plot_timecourse_eeg_emg(config{ipatient}, dat_LFP, ipart, imarker, 'eeg',config{ipatient}.align.channel{imarker},config{ipatient}.epoch.toi{imarker}, true);
%                     
%                     dtx_plot_timecourse_eeg_emg(config{ipatient}, dat_LFP, ipart, imarker, 'eeg',config{ipatient}.LFP.motorcortex{imarker},config{ipatient}.epoch.toi{imarker}, true);
%                     
%                     dtx_plot_overdraw_allchannels(config{ipatient},dat_LFP, ipart, imarker, true);
%                     
%                     dtx_plot_avg_allchannels(config{ipatient},dat_LFP,ipart,imarker,true);
%                     
%                     dtx_plot_SlowWaveTopographyTimecouse(config{ipatient},dat_LFP, ipart, imarker, true);
%                     
                    %[TFR_eeg{ipatient}] = dtx_plot_TFR(config{ipatient}, dat_LFP, ipart, imarker, true);
                    
                    %keep eeg for plot of both SW_R and SW_L in the same figure
                    iEEG = iEEG + 1;
                    dataEEG{iEEG} = dat_LFP{ipart}{imarker};
                    dataEEG{iEEG}.LFP.name = config{ipatient}.LFP.name{imarker};
                    
                    %get patient data for figure pool
                    if ipart == length(config{ipatient}.directorylist)
                        [data_avg_allchans{imarker}{ipatient},...
                            data_avg_chanalign{imarker}{ipatient},...
                            data_avg_EMG{imarker}{ipatient}] =...
                            dtx_get_patient_data(config{ipatient}, dat_LFP, ipart, imarker);
                    end
                    
%                 end %if is EMG
            else %if dat_LFP{ipart}{imarker} is empty
                data_avg_allchans{imarker}{ipatient} = [];
                data_avg_chanalign{imarker}{ipatient} = [];
                data_avg_EMG{imarker}{ipatient} = [];                
            end
            
        end %imarker
        
        %1 figure common for all markers (SW_R and SW_L)
        
%         dtx_plot_SlowWaveTopography(config{ipatient},dataEEG,ipart, true);
        
        if config{ipatient}.continous == true
            Seizure_Infos{ipatient} =...
                dtx_plot_count_seizure(config{ipatient}, MuseStruct, ipart, false);
        else
            Seizure_Infos{ipatient} = [];
        end
        
    end %ipart
    clear dat_LFP dataEEG MuseStruct
end %ipatient


%% save data all patients

%save. Same datasavedir for all patients.
%data are the data of the last ipart : so the 'merged' data
save(fullfile(config{1}.datasavedir,'All_patients_data_avg_allchansn.mat'),'data_avg_allchans');
save(fullfile(config{1}.datasavedir,'All_patients_data_avg_chanalign.mat'),'data_avg_chanalign');
save(fullfile(config{1}.datasavedir,'All_patients_data_avg_EMG.mat'),'data_avg_EMG');
save(fullfile(config{1}.datasavedir,'All_patients_Seizure_Infos.mat'),'Seizure_Infos');


%% Figures all patients
for do_normalize = [true false]
    
    [config] = dtx_setparams_patients_lgi1([]);
    
    % Same datasavedir for all patients
    load(fullfile(config{1}.datasavedir,'All_patients_data_avg_allchansn.mat'),'data_avg_allchans');
    load(fullfile(config{1}.datasavedir,'All_patients_data_avg_chanalign.mat'),'data_avg_chanalign');
    load(fullfile(config{1}.datasavedir,'All_patients_data_avg_EMG.mat'),'data_avg_EMG');
    load(fullfile(config{1}.datasavedir,'All_patients_Seizure_Infos.mat'),'Seizure_Infos');
    
    % keep same data structure as for patient by patient analysis
    ipart = 1;
    ipatient = 1;
    
    %set config
    config{ipatient}.prefix = 'AllPatients-';
    config{ipatient}.imagesavedir = [config{ipatient}.imagesavedir,'\..\AllPatients'];
    config{ipatient}.LFP.emg                   = {'no','no','EMG','EMG'};%same index as associated EEG. 'no' if no EMG associated to this seizure side
    config{ipatient}.merge = false;
    
    
    % Append data : one patient is one trial
    for imarker = 1:length(config{ipatient}.LFP.name)
        
        isdata = [];
        isdata = ~cellfun('isempty',data_avg_allchans{imarker}); %remove empty arrays for append
        if any(isdata)
            data_avg_allchans_allpatients{ipart}{imarker}                = ft_appenddata([],data_avg_allchans{imarker}{isdata});
        else
            data_avg_allchans_allpatients{ipart}{imarker}                = [];
        end
        
        isdata = [];
        isdata = ~cellfun('isempty',data_avg_chanalign{imarker}); %remove empty arrays for append
        if any(isdata)
            data_avg_chanalign_allpatients{ipart}{imarker}                  = ft_appenddata([],data_avg_chanalign{imarker}{isdata});
        else
            data_avg_chanalign_allpatients{ipart}{imarker}                = [];
        end
        
        isdata = [];
        isdata = ~cellfun('isempty',data_avg_EMG{imarker}); %remove empty arrays for append
        if any(isdata)
            data_avg_EMG_allpatients{ipart}{imarker}        = ft_appenddata([],data_avg_EMG{imarker}{isdata});
        else
            data_avg_EMG_allpatients{ipart}{imarker}                = [];
        end
        
    end
    
    
    %normalize by dividing by value of max SlowWave at t=0 (alignment).
    %normalize emg by dividing by max(EMG) beetween 0 and 2
    if do_normalize
        config{ipatient}.prefix = 'AllPatients-normalized-';
        config{ipatient}.imagesavedir = [config{ipatient}.imagesavedir,'-normalized'];
        
        for imarker = 1:length(config{ipatient}.LFP.name)
            
            for itrial = 1:length(config) %one trial is one patient
                %set norm values
                %norm_eeg : value at t=0 for align channel for eeg, or for
                %c4c3 channel for emg
                %norm_emg : max of trial between 0 and 2s
                t = [];
                norm_eeg = 0;
                norm_emg = 0;
                if ~isempty(data_avg_chanalign_allpatients{ipart}{imarker})
                    if itrial <= length(data_avg_chanalign_allpatients{ipart}{imarker}.trial)
                        t = data_avg_chanalign_allpatients{ipart}{imarker}.time{itrial};
                        norm_eeg = data_avg_chanalign_allpatients{ipart}{imarker}.trial{itrial}(t==0);
                    end
                end
                if ~isempty(data_avg_EMG_allpatients{ipart}{imarker})
                    if itrial <= length(data_avg_EMG_allpatients{ipart}{imarker}.trial)
                        motorcortex_indx = find(strcmp(data_avg_EMG_allpatients{ipart}{imarker}.label, config{ipatient}.LFP.motorcortex{imarker}));
                        t = data_avg_EMG_allpatients{ipart}{imarker}.time{itrial};
                        norm_eeg = data_avg_EMG_allpatients{ipart}{imarker}.trial{itrial}(motorcortex_indx,t==0);
                        norm_emg = max(data_avg_EMG_allpatients{ipart}{imarker}.trial{itrial}(end,t>=0 & t<=2)); %end = EMG, defined by dtx_get_patient_data
                    end
                end
                
                %apply norm values to data
                if ~isempty(data_avg_allchans_allpatients{ipart}{imarker})
                    if itrial <= length(data_avg_allchans_allpatients{ipart}{imarker}.trial)
                        data_avg_allchans_allpatients{ipart}{imarker}.trial{itrial} = data_avg_allchans_allpatients{ipart}{imarker}.trial{itrial} / norm_eeg;
                    end
                end
                
                if ~isempty(data_avg_EMG_allpatients{ipart}{imarker})
                    if itrial <= length(data_avg_EMG_allpatients{ipart}{imarker}.trial)
                        data_avg_EMG_allpatients{ipart}{imarker}.trial{itrial}(motorcortex_indx,:) = data_avg_EMG_allpatients{ipart}{imarker}.trial{itrial}(motorcortex_indx,:) / norm_eeg;
                        data_avg_EMG_allpatients{ipart}{imarker}.trial{itrial}(end,:) = data_avg_EMG_allpatients{ipart}{imarker}.trial{itrial}(end,:) / norm_emg;
                    end
                end
                
                if ~isempty(data_avg_chanalign_allpatients{ipart}{imarker})
                    if itrial <= length(data_avg_chanalign_allpatients{ipart}{imarker}.trial)
                        data_avg_chanalign_allpatients{ipart}{imarker}.trial{itrial} = data_avg_chanalign_allpatients{ipart}{imarker}.trial{itrial} / norm_eeg;
                    end
                end
                
            end
        end %imarker
        
    end %do_normalize
    
    
    %plot data : exactly the same protocol as patient by patient
    dataEEG = []; %for figure with right and left eeg data, without emg data
    iEEG = 0;
    
    for imarker = 1:length(config{ipatient}.LFP.name)
        
        %EMG analysis
        if ~isempty(data_avg_EMG_allpatients{ipart}{imarker})
            
            emgToPlot = 'EMG';
            
            dtx_plot_emg_method(config{ipatient},data_avg_EMG_allpatients,ipart,imarker,emgToPlot,true);
            
            dtx_plot_timecourse_eeg_emg(config{ipatient}, data_avg_EMG_allpatients, ipart, imarker, 'emg',emgToPlot,config{ipatient}.epoch.toi{imarker}, true);
            
            dtx_plot_comparison_eeg_emg(config{ipatient}, data_avg_EMG_allpatients, ipart, imarker, true, config{ipatient}.LFP.motorcortex{imarker}, emgToPlot);
            
            dtx_plot_comparison_eeg_emg_trialbytrial(config{ipatient}, data_avg_EMG_allpatients, ipart, imarker, true, config{ipatient}.LFP.motorcortex{imarker}, emgToPlot)
            
        end
        
        %EEG analysis
        
        if ~isempty(data_avg_allchans_allpatients{ipart}{imarker})
                        
            dtx_plot_timecourse_eeg_emg(config{ipatient}, data_avg_chanalign_allpatients, ipart, imarker, 'eeg','chan_align',config{ipatient}.epoch.toi{imarker}, true);
            
            dtx_plot_timecourse_eeg_emg(config{ipatient}, data_avg_allchans_allpatients, ipart, imarker, 'eeg',config{ipatient}.LFP.motorcortex{imarker},config{ipatient}.epoch.toi{imarker}, true);
            
            dtx_plot_overdraw_allchannels(config{ipatient},data_avg_allchans_allpatients, ipart, imarker, true);
            
            dtx_plot_avg_allchannels(config{ipatient},data_avg_allchans_allpatients,ipart,imarker,true);
            
            dtx_plot_SlowWaveTopographyTimecouse(config{ipatient},data_avg_allchans_allpatients, ipart, imarker, true);
            
            iEEG = iEEG + 1;
            dataEEG{iEEG} = data_avg_allchans_allpatients{ipart}{imarker};
            dataEEG{iEEG}.LFP.name = config{ipatient}.LFP.name{imarker};
            
        end
        
    end %imarker
    
    %1 figure common for all markers (SW_R and SW_L)
    dtx_plot_SlowWaveTopography(config{ipatient},dataEEG,ipart, true);
    
end %do_normalize
%AJOUTER CV CV2 FANO FACTOR DE CHAQUE PATIENT SUR UN PLOT



