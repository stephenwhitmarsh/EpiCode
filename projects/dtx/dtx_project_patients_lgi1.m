
%% Analyse of LGI1 patients seizures



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
            
            
            [MuseStruct]    = readMuseMarkers(config{ianalyse}, false);
            [MuseStruct]    = alignMuseMarkers(config{ianalyse},MuseStruct, false);
            [dat_LFP]       = readLFP(config{ianalyse}, MuseStruct, false, false); %dat_LFP{ipart}{imarker}
            
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
        
        
        %% Choose EEG and EMG channels according to the setparams
        iEEG = [];
        iEMG = [];
        
        for ichannel = 1 : length(dat_LFP{ipart}{imarker}.label)
            for imarker = 1:length(config{ipatient}.LFP.name)
                if strcmp(dat_LFP{ipart}{imarker}.label{ichannel},config{ipatient}.align.channel{imarker})
                    iEEG{imarker} = ichannel;
                end
                if strcmp(dat_LFP{ipart}{imarker}.label{ichannel},config{ipatient}.LFP.emg{imarker})
                    iEMG{imarker} = ichannel;
                end
                if strcmp(config{ipatient}.LFP.emg{imarker},'no')
                    iEMG{imarker} = {'no'};
                end
            end
            
        end
        if length(iEMG)==1
            iEMG{2}={'no'};
        end
        %% Plot
        
        fprintf('***********************************\n');
        fprintf('***********************************\n');
        fprintf('**** Plot data, for patient %d ****\n',ipatient);
        fprintf('***********************************\n');
        fprintf('***********************************\n\n');
        
        
        %2 figures : For SW_R and SW_L
        for imarker = 1:length(config{ipatient}.LFP.name)
            
            dtx_plot_timecourse_eeg_emg(config{ipatient}, dat_LFP{ipart}, iEEG, iEMG, imarker,true);
            
            dtx_plot_TFR(config{ipatient}, dat_LFP{ipart}, iEEG, imarker, true);
            
            if ~any(strcmp(iEMG{imarker},'no') || iEMG{imarker} == true)
                
                dtx_plot_comparison_eeg_emg(config{ipatient}, dat_LFP{ipart}, iEEG, iEMG, imarker,true); %find which strategy of EMG quantif
                
            end
            
        end
        
        %1 figure common for SW_R and SW_L
        dtx_plot_overdraw_allchannels(config{ipatient},dat_LFP{ipart},iEEG,iEMG,true,true);
        
        dtx_plot_SlowWaveTopography(config{ipatient},dat_LFP{ipart},true,true);
        
        dtx_plot_SlowWaveTopographyTimecouse(config{ipatient},dat_LFP{ipart},true);
        
        dtx_plot_count_seizure(config{ipatient}, MuseStruct, true)
        
    end %ipatient
end %merge_eeg


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% average over trials for plotting of all channels
% for abscisse_max = 15%[2, 5, 10]
%
%     cfgtemp                 = [];
%     cfgtemp.vartrllength    = 2;
%     %dat_micro_rptavg        = ft_timelockanalysis(cfgtemp,dat_micro{1}{1});
%     dat_LFP_rptavg        = ft_timelockanalysis(cfgtemp,dat_LFP{1}{1});
%
%     fig=figure;
%     hold;
%     h=200;
%     for i = 1 : size(dat_LFP_rptavg.label,1)
%         %subplot(size(dat_LFP_rptavg.label,1),1,size(dat_LFP_rptavg.label,1)-i+1);
%         %         plot(dat_micro_rptavg.time,dat_micro_rptavg.avg)
%         plot(dat_LFP_rptavg.time,dat_LFP_rptavg.avg(i,:)+h*i,'color','black');
%
%     end
%     title(sprintf('Average of all SlowWaves of %s',config{ipatient}.prefix(1:end-1)),'Interpreter','none');
%     axis tight;
%     xlabel('Time (s)');
%     ylabel('Channel name');
%     tick = h;
%     yticks(0 : tick : length(dat_LFP_rptavg.label)*h);
%     set(gca, 'YTickLabel',[dat_LFP_rptavg.label(end)', dat_LFP_rptavg.label(1:end-1)']); %why Y axis not in the good order
%     set(gca,'TickDir','out');
%     %         if abscisse_max == 2
%     %             xlim([-2, 2]);
%     %         else
%     %             xlim([-5, abscisse_max]);
%     %         end
%     %         line ([abscisse_max-0.5,abscisse_max-0.5],[h*i-2500,h*i-500],'color','r','LineWidth',2);
%     %         text(abscisse_max-0.45, h*i-1400,'2mV','FontSize',15);
%
%
%     %         % print to file
%     %         fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
%     %         set(fig,'PaperOrientation','landscape');
%     %         set(fig,'PaperUnits','normalized');
%     %         set(fig,'PaperPosition', [0 0 1 1]);
%     %         print(fig, '-dpdf', fullfile(config{ipatient}.imagesavedir,'average_SlowWave',[config{ipatient}.prefix,'_Avg_SlowWaves_MACRO_',num2str(abscisse_max),'s']),'-r600');
%     %         print(fig, '-dpng', fullfile(config{ipatient}.imagesavedir,'average_SlowWave',[config{ipatient}.prefix,'_Avg_SlowWaves_MACRO_',num2str(abscisse_max),'s']),'-r600');
%     %
%     %         close all
%
%
% end
