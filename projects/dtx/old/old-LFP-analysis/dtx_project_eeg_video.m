%% Analyse of LGI1 patients seizures

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
feature('DefaultCharacterSet', 'CP1252') % To fix bug for weird character problems in reading neuralynx


config = dtx_setparams_eegvideo([]);
%cfg=config{1};


for irat = 1:length(config)
    
    for itemp = 1:3
        if itemp == 1
            config = dtx_setparams_eegvideo([]);
            config{1}.imagesavedir = fullfile(config{1}.imagesavedir, 'hp_0.1');
            config{1}.LFP.hpfilter = 'yes';
            config{1}.LFP.hpfreq = 0.1;
        end
        if itemp == 2
            config = dtx_setparams_eegvideo([]);
            config{1}.imagesavedir = fullfile(config{1}.imagesavedir, 'hp_0.5');
            config{1}.LFP.hpfilter = 'yes';
            config{1}.LFP.hpfreq = 0.5;
        end
        if itemp == 3
            config = dtx_setparams_eegvideo([]);
            config{1}.imagesavedir = fullfile(config{1}.imagesavedir, 'nohp');
            config{1}.LFP.hpfilter = 'no';
        end
        
    %% Get LFP data
    %config{1}.imagesavedir = fullfile(config{1}.imagesavedir, 'hp_1');
    [MuseStruct]               = readMuseMarkers(config{irat}, false);
    
    %check_nr_crises_startend(config{irat}, MuseStruct, 1);
    
    [MuseStruct]               = alignMuseMarkers(config{irat},MuseStruct, false);
    [dat_LFP]                  = readLFP(config{irat}, MuseStruct, false, false);
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
    
    [datatest, MuseStructtest] = removetrial_marker(config{irat}, dat_LFP, MuseStruct,'keep','SlowWave_EMG','SlowWave_EMG',0,0,-5,5, true);

    [datatest, MuseStructtest] = removetrial_marker(config{irat}, dat_LFP, MuseStruct,'remove','BAD__START__','BAD__END__',0,0,0,0, true);
    [dat_LFP, MuseStruct] = dtx_removeartefactsLFP(dat_LFP, MuseStruct, config{irat}, true);
    
    
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
                
                %is EMG
                if any(strfind(config{irat}.LFP.name{imarker}, 'EMG'))
                    
                    emgToPlot = config{irat}.LFP.emg{imarker};
                    
                    dtx_plot_emg_method(config{irat},dat_LFP,ipart,imarker,emgToPlot,true);
                    
                    dtx_plot_timecourse_eeg_emg(config{irat}, dat_LFP, ipart, imarker, 'emg',emgToPlot,config{irat}.EMG.toi, true);
                    
                    dtx_plot_comparison_eeg_emg(config{irat}, dat_LFP, ipart, imarker, true, config{irat}.align.channel{imarker},emgToPlot);
                    
                    %is EEG
                else
                    
                    dtx_plot_timecourse_eeg_emg(config{irat}, dat_LFP, ipart, imarker, 'eeg',config{irat}.align.channel{imarker},config{irat}.epoch.toi{imarker}, true);
                    
                    dtx_plot_overdraw_allchannels(config{irat},dat_LFP, ipart, imarker, true);
                    
                    dtx_plot_avg_allchannels(config{irat},dat_LFP,ipart,imarker,true);
                    
                    [TFR_eeg{irat}] = dtx_plot_TFR(config{irat}, dat_LFP, ipart, imarker, true);
                    
                    %keep eeg for plot of both SW_R and SW_L in the same figure
                    iEEG = iEEG + 1;
                    dataEEG{iEEG} = dat_LFP{ipart}{imarker};
                    dataEEG{iEEG}.LFP.name = config{irat}.LFP.name{imarker};
                    
                    %get rat data for figure pool
                    [data_avg_allchans{imarker}{irat},...
                        data_avg_chanalign{imarker}{irat},...
                        data_avg_EMG{imarker}{irat}] =...
                        dtx_get_patient_data(config{irat}, dat_LFP, ipart, imarker);
                    
                end %                 end %if is EMG
            else %if dat_LFP{ipart}{imarker} is empty
                data_avg_allchans{imarker}{irat} = [];
                data_avg_chanalign{imarker}{irat} = [];
                data_avg_EMG{imarker}{irat} = [];
            end
            
        end %imarker
        
        %1 figure common for all markers (SW_R and SW_L)
        
        if config{irat}.continous == true
            Seizure_Infos{irat} =...
                dtx_plot_count_seizure(config{irat}, MuseStruct, ipart, true);
        else
            Seizure_Infos{irat} = [];
        end
        
    end %ipart
    clear dat_LFP dataEEG MuseStruct
    end
end %irat

%% save data all rats

%save. Same datasavedir for all patients.
%data are the data of the last ipart : so the 'merged' data
save(fullfile(config{1}.datasavedir,'All_rats_data_avg_allchansn.mat'),'data_avg_allchans');
save(fullfile(config{1}.datasavedir,'All_rats_data_avg_chanalign.mat'),'data_avg_chanalign');
save(fullfile(config{1}.datasavedir,'All_rats_data_avg_EMG.mat'),'data_avg_EMG');
save(fullfile(config{1}.datasavedir,'All_rats_Seizure_Infos.mat'),'Seizure_Infos');
save(fullfile(config{1}.datasavedir,'All_rats_TFR_eeg.mat'),'TFR_eeg');