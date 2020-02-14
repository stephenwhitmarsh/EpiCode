
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

for merge_eeg = [false true] %one analysis patient per patient, one other merging
    

    if merge_eeg == false
        nb_of_patients = length(config);
    elseif merge_eeg == true
        nb_of_patients = length(mergeindex);
    end
    
    
    for ipatient = 1:nb_of_patients %go patient/eeg by patient/eg

        
        %% Get LFP data
        
        if merge_eeg == false
            idata = ipatient;
        elseif merge_eeg == true
            idata = mergeindex{ipatient};
        end
        
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
                data_temp{ianalyse} = dat_LFP;
                MuseStruct_temp{ianalyse} = MuseStruct;
            end
        end
        
        if merge_eeg == true 
            [config{ipatient}, MuseStruct, dat_LFP_EMG, dat_LFP] = dtx_merge_patientslgi1(mergeindex,ipatient,MuseStruct_temp,data_temp);
        end
        

        %% Plot
        
        fprintf('***********************************\n');
        fprintf('***********************************\n');
        fprintf('**** Plot data, for %s ****\n',config{ipatient}.prefix(1:end-1));
        fprintf('***********************************\n');
        fprintf('***********************************\n\n');
        
        
        %1 figure per marker type
        for imarker = 1:length(config{ipatient}.LFP.name)
            
            dtx_plot_timecourse_eeg_emg(config{ipatient}, dat_LFP{ipart},imarker, 'eeg', true);
            dtx_plot_timecourse_eeg_emg(config{ipatient}, dat_LFP_EMG{ipart},imarker, 'emg', true);
            
            dtx_plot_TFR(config{ipatient}, dat_LFP{ipart}, imarker, true);
            
            dtx_plot_overdraw_allchannels(config{ipatient},dat_LFP{ipart},imarker,true);
            
            dtx_plot_avg_allchannels(config{ipatient},dat_LFP{ipart},imarker,true);
            
            dtx_plot_SlowWaveTopographyTimecouse(config{ipatient},dat_LFP{ipart}, imarker, true);
            
            if isfield(config{ipatient}.LFP,'emg')
                if ~(strcmp(config{ipatient}.LFP.emg{imarker},'no'))
                    if ~(config{ipatient}.LFP.emg{imarker} == false)
                        
                        dtx_plot_comparison_eeg_emg(config{ipatient}, dat_LFP_EMG{ipart}, imarker,true); %find which strategy of EMG quantif
                    
                    end
                end
            end
            
        end
        
        %1 figure common for all markers (SW_R and SW_L)
        
        dtx_plot_SlowWaveTopography(config{ipatient},dat_LFP{ipart},true,true);
        
        dtx_plot_count_seizure(config{ipatient}, MuseStruct, true)
        
    end %ipatient
end %merge_eeg

