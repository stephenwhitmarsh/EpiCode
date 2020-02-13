
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
merge_eeg = false;


for ipatient = 1:length(config)
    
    
    %% Get LFP data
    
    
    [MuseStruct]    = readMuseMarkers(config{ipatient}, false);
    [MuseStruct]    = alignMuseMarkers(config{ipatient},MuseStruct, false);
    [dat_LFP]       = readLFP(config{ipatient}, MuseStruct, false, false); %dat_LFP{ipart}{imarker}
    
    %inverse the data according to standard clinical visualisation
    for imarker = 1:size(dat_LFP{1},2)
        for itrial = 1:size(dat_LFP{1}{imarker}.trial,2)
            dat_LFP{1}{imarker}.trial{itrial} = -dat_LFP{1}{imarker}.trial{itrial};
        end
    end
   
    [config{ipatient}, dat_LFP] = dtx_merge_patientslgi1(ipatient_to_merge,dat_LFP);
    
    %% Choose EEG and EMG channels according to the setparams
    for ichannel = 1 : length(dat_LFP{1}{1}.label)
        for imarker = 1:length(config{ipatient}.LFP.name)
            
            if strcmp(dat_LFP{1}{imarker}.label{ichannel},config{ipatient}.align.channel{imarker})
                iEEG{imarker} = ichannel;
            end
            if strcmp(dat_LFP{1}{imarker}.label{ichannel},config{ipatient}.LFP.emg{imarker})
                iEMG{imarker} = ichannel;
            end
            if strcmp(config{ipatient}.LFP.emg{imarker},'no')
                iEMG{imarker} = {'no'};
            end
        end
    end
    
    %% Plot    
    
    fprintf('***********************************\n');
    fprintf('***********************************\n');
    fprintf('**** Plot data, for patient %d ****\n',ipatient);
    fprintf('***********************************\n');
    fprintf('***********************************\n\n');
    
    
    %2 figures : For SW_R and SW_L
    for imarker = 1:length(config{ipatient}.LFP.name)
        
        dtx_plot_timecourse_eeg_emg(config{ipatient}, dat_LFP{1}, iEEG, iEMG, imarker,false);
        
        dtx_plot_TFR(config{ipatient}, dat_LFP{1}, iEEG, imarker, false);
        
        if ~any(strcmp(iEMG{imarker},'no') || iEMG{imarker} == false)
            
            dtx_plot_comparison_eeg_emg(config{ipatient}, dat_LFP{1}, iEEG, iEMG, imarker,true); %find which strategy of EMG quantif
            
        end
        
    end
    
    %1 figure common for SW_R and SW_L
    dtx_plot_overdraw_allchannels(config{ipatient},dat_LFP{1},iEEG,iEMG,true,true);
    
    dtx_plot_SlowWaveTopography(config{ipatient},dat_LFP{1},true,true);
    
    dtx_plot_SlowWaveTopographyTimecouse(config{ipatient},dat_LFP{1},true);
    
    dtx_plot_count_seizure(config{ipatient}, MuseStruct, true)
    
    
end


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
