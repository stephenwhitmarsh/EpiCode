
%% Analyse manips WOD
%Paul Baudin

%% Set parameters
if ispc
    addpath \\lexport\iss01.charpier\analyses\lgi1\Git-Paul\fieldtrip;
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\projects\wod'));
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\external'));
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\shared'));
    addpath \\lexport\iss01.charpier\analyses\lgi1\Git-Paul\CEDMATLAB\CEDS64ML
    CEDS64LoadLib('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\CEDMATLAB\CEDS64ML');

elseif isunix
    addpath /network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/fieldtrip;
    addpath (genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/projects/wod'));
    addpath (genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/external'));
    addpath (genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/shared'));
end
ft_defaults
feature('DefaultCharacterSet', 'CP1252') % To fix bug for weird character problems in reading neuralynx

[config] = wod_setparams([]);
cfg=config{1};

for irat = 1:length(config)
    ipart = 1;
    
    %% Get LFP data
    
% 	For readCEDMarkers and CED2mat : Can only be done on windows because 
%   of CEDS64 library. For Linux, set force argument to false, it will 
%   load the results.
    [MuseStruct]               = readCEDMarkers(config{irat}, false);
    CEDcontinuous_to_mat(config{irat},false);
    [MuseStruct]               = alignMuseMarkers(config{irat},MuseStruct, true);
    
    
    
%     %write markers to Muse Neuralynx to see if they are coherents
%     for idir = 1:length(MuseStruct{ipart})
%         writeMuseMarkerfile(MuseStruct,fullfile(config{irat}.rawdir,config{irat.directorylist{idir},'mrkMusetest.mrk'));
%     end
    
%créer un marker avec toutes les WOD (de irat ipart idir) avec comme info
%nom marker initial, nom canal associé
%marker.synctime
%marker.origname


    %aligner chaque WOD à son pic avec le marker_total créé.
    [MuseStruct]               = alignMuseMarkers(config{irat},MuseStruct, false);
    
    %read avec chaque marker séparé (défini dans muse.startend
    [dat_LFP]                  = readLFP(config{irat}, MuseStruct, false, false);
    %dat_LFP{ipart}{imarker}
       
    %append data de chaque marker, 
    
    [dat_LFP, MuseStruct] = dtx_removeartefactsLFP(dat_LFP, MuseStruct, config{irat}, false);
    
    
    %% Plot
    
    
    
    for ipart = 1:length(dat_LFP)
        
        fprintf('*********************************************\n');
        fprintf('*********************************************\n');
        fprintf('**** Plot data, for %s, part %d ****\n',config{irat}.prefix(1:end-1), ipart);
        fprintf('*********************************************\n');
        fprintf('*********************************************\n\n');
        
        
        dataEEG = []; %for figure with right and left eeg data, without emg data
        iEEG = 0;
        
        %1 figure per marker type
        for imarker = 1:length(config{irat}.LFP.name)

            if ~isempty(dat_LFP{ipart}{imarker})
                
         
                                        
                    dtx_plot_timecourse_eeg_emg(config{irat}, dat_LFP, ipart, imarker, 'eeg',config{irat}.align.channel{imarker},config{irat}.epoch.toi{imarker}, true);
                    
                    dtx_plot_timecourse_eeg_emg(config{irat}, dat_LFP, ipart, imarker, 'eeg',config{irat}.LFP.motorcortex{imarker},config{irat}.epoch.toi{imarker}, true);
                    
                    dtx_plot_overdraw_allchannels(config{irat},dat_LFP, ipart, imarker, true);
                    
                    dtx_plot_avg_allchannels(config{irat},dat_LFP,ipart,imarker,true);
                    
                    dtx_plot_SlowWaveTopographyTimecouse(config{irat},dat_LFP, ipart, imarker, true);
                     
                    %[TFR_eeg{ipatient}] = dtx_plot_TFR(config{ipatient}, dat_LFP, ipart, imarker, true);
                    
                    %keep eeg for plot of both SW_R and SW_L in the same figure
                    iEEG = iEEG + 1;
                    dataEEG{iEEG} = dat_LFP{ipart}{imarker};
                    dataEEG{iEEG}.LFP.name = config{irat}.LFP.name{imarker};
                    
                    %get patient data for figure pool
                    if ipart == length(config{irat}.directorylist)
                        [data_avg_allchans{imarker}{irat},...
                            data_avg_chanalign{imarker}{irat},...
                            data_avg_EMG{imarker}{irat}] =...
                            dtx_get_patient_data(config{irat}, dat_LFP, ipart, imarker);
                    end
                    
%                 end %if is EMG
            else %if dat_LFP{ipart}{imarker} is empty
                data_avg_allchans{imarker}{irat} = [];
                data_avg_chanalign{imarker}{irat} = [];
                data_avg_EMG{imarker}{irat} = [];                
            end
            
        end %imarker
        
        %1 figure common for all markers (SW_R and SW_L)
        
%         dtx_plot_SlowWaveTopography(config{ipatient},dataEEG,ipart, true);
        
        if config{irat}.continous == true
            Seizure_Infos{irat} =...
                dtx_plot_count_seizure(config{irat}, MuseStruct, ipart, false);
        else
            Seizure_Infos{irat} = [];
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
    irat = 1;
    
    %set config
    config{irat}.prefix = 'AllPatients-';
    config{irat}.imagesavedir = [config{irat}.imagesavedir,'\..\AllPatients'];
    config{irat}.LFP.emg                   = {'no','no','EMG','EMG'};%same index as associated EEG. 'no' if no EMG associated to this seizure side
    config{irat}.merge = false;
    
    
    % Append data : one patient is one trial
    for imarker = 1:length(config{irat}.LFP.name)
        
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
        config{irat}.prefix = 'AllPatients-normalized-';
        config{irat}.imagesavedir = [config{irat}.imagesavedir,'-normalized'];
        
        for imarker = 1:length(config{irat}.LFP.name)
            
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
                        motorcortex_indx = find(strcmp(data_avg_EMG_allpatients{ipart}{imarker}.label, config{irat}.LFP.motorcortex{imarker}));
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
    
    for imarker = 1:length(config{irat}.LFP.name)
        
        %EMG analysis
        if ~isempty(data_avg_EMG_allpatients{ipart}{imarker})
            
            emgToPlot = 'EMG';
            
            dtx_plot_emg_method(config{irat},data_avg_EMG_allpatients,ipart,imarker,emgToPlot,true);
            
            dtx_plot_timecourse_eeg_emg(config{irat}, data_avg_EMG_allpatients, ipart, imarker, 'emg',emgToPlot,config{irat}.epoch.toi{imarker}, true);
            
            dtx_plot_comparison_eeg_emg(config{irat}, data_avg_EMG_allpatients, ipart, imarker, true, config{irat}.LFP.motorcortex{imarker}, emgToPlot);
            
            dtx_plot_comparison_eeg_emg_trialbytrial(config{irat}, data_avg_EMG_allpatients, ipart, imarker, true, config{irat}.LFP.motorcortex{imarker}, emgToPlot)
            
        end
        
        %EEG analysis
        
        if ~isempty(data_avg_allchans_allpatients{ipart}{imarker})
                        
            dtx_plot_timecourse_eeg_emg(config{irat}, data_avg_chanalign_allpatients, ipart, imarker, 'eeg','chan_align',config{irat}.epoch.toi{imarker}, true);
            
            dtx_plot_timecourse_eeg_emg(config{irat}, data_avg_allchans_allpatients, ipart, imarker, 'eeg',config{irat}.LFP.motorcortex{imarker},config{irat}.epoch.toi{imarker}, true);
            
            dtx_plot_overdraw_allchannels(config{irat},data_avg_allchans_allpatients, ipart, imarker, true);
            
            dtx_plot_avg_allchannels(config{irat},data_avg_allchans_allpatients,ipart,imarker,true);
            
            dtx_plot_SlowWaveTopographyTimecouse(config{irat},data_avg_allchans_allpatients, ipart, imarker, true);
            
            iEEG = iEEG + 1;
            dataEEG{iEEG} = data_avg_allchans_allpatients{ipart}{imarker};
            dataEEG{iEEG}.LFP.name = config{irat}.LFP.name{imarker};
            
        end
        
    end %imarker
    
    %1 figure common for all markers (SW_R and SW_L)
    dtx_plot_SlowWaveTopography(config{irat},dataEEG,ipart, true);
    
end %do_normalize
%AJOUTER CV CV2 FANO FACTOR DE CHAQUE PATIENT SUR UN PLOT



