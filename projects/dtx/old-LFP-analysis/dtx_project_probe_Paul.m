
%% Analyse of DTX seizures


% Stephen Whitmarsh, modified by Paul Baudin
%
% requires bandpassFilter.m from Mario
% requires releaseDec2015 from Neuralynx website
%requires signal processing toolbox


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

config = dtx_setparams_probe([]);

%% LFP analysis
% 
% for irat = 1:length(config)
%     %% Get right LFP data
%     % read muse markers
%     
%     [MuseStruct]    = readMuseMarkers(config{irat}, true);
%     
%     % align Muse markers according to peaks and detect whether they contain artefacts
%     [MuseStruct]    = alignMuseMarkers(config{irat},MuseStruct, true);
%     
%     %[MuseStruct] =   dtx_remove_wrong_seizure(config{irat},
%     %MuseStruct,false, false); %only for spike analysis
%     [dat_LFP] = readLFP(config{irat}, MuseStruct, true, true);
%     
%     [dat_LFP, MuseStruct] = dtx_removeartefactsLFP(dat_LFP, MuseStruct, config{irat}, true);
%     %dat_LFP{ipart}{imarker}
%     
%     
%     ipart = 1;
%     
%     for imarker = 1:size(dat_LFP{ipart},2)
%         
%         if irat == 2 %wrong channel name during acquisition
%             
%             if strcmp(config{irat}.LFP.electrodetoplot{imarker},'ECoGS1')
%                 config{irat}.LFP.electrodetoplot{imarker} = 'ECoGM1G';
%             end
%             
%             for ilabel = 1:length(dat_LFP{ipart}{imarker}.label)
%                 if strcmp(dat_LFP{ipart}{imarker}.label{ilabel}, 'ECoGM1')
%                     dat_LFP{ipart}{imarker}.label{ilabel} = 'ECoGM1D';
%                 end
%                 if strcmp(dat_LFP{ipart}{imarker}.label{ilabel}, 'ECoGS1')
%                     dat_LFP{ipart}{imarker}.label{ilabel} = 'ECoGM1G';
%                 end
%                 if strcmp(config{irat}.labels.macro{ilabel}, 'ECoGM1')
%                     config{irat}.labels.macro{ilabel} = 'ECoGM1D';
%                 end
%                 if strcmp(config{irat}.labels.macro{ilabel}, 'ECoGS1')
%                     config{irat}.labels.macro{ilabel} = 'ECoGM1G';
%                 end
%             end
%             
%         end
%         
%         electrodeToPlot = config{irat}.LFP.electrodetoplot{imarker};
%         
%         dtx_plot_timecourse_eeg_emg(config{irat}, dat_LFP, ipart, imarker, 'eeg',electrodeToPlot,[-inf, inf], true);
%         dtx_plot_timecourse_eeg_emg(config{irat}, dat_LFP, ipart, imarker, 'eeg',electrodeToPlot,[-5, 25], true);
%         dtx_plot_timecourse_eeg_emg(config{irat}, dat_LFP, ipart, imarker, 'eeg',electrodeToPlot,[-5, 5], true);
%         dtx_plot_TFR(config{irat}, dat_LFP, ipart, imarker, false, true);
%         
%         [data_avg_allchans{imarker}{irat}, data_avg_chantoplot{imarker}{irat}, data_TFR{imarker}{irat}, ~] =...
%             dtx_get_patient_data(config{irat}, dat_LFP, ipart, imarker);
%         
%         Seizure_Infos{imarker}{irat} = dtx_plot_count_seizure(config{irat}, MuseStruct, ipart, true); %need to be done with all markers, not markers removed wrong seizures
%         
%     end
%     
%     %figure common for 2 LFP names
%     imarker = 1;
%     dtx_plot_avg_allchannels(config{irat},dat_LFP,ipart,imarker,config{irat}.epoch.toi{imarker},true);
%     dtx_plot_avg_allchannels(config{irat},dat_LFP,ipart,imarker,[-2 2],true);
%     dtx_plot_overdraw_allchannels(config{irat},dat_LFP, ipart, imarker, true);
%     dtx_plot_comparison_severaleeg(config{irat},dat_LFP,ipart,imarker,{config{irat}.LFP.electrodetoplot{:}}, [-2 2], true, true)
%     
% end
% 
% %% save data all rats
% 
% %save. Same datasavedir for all rats.
% %data are the data of the last ipart : so the 'merged' data
% save(fullfile(config{1}.datasavedir,'All_rats_data_avg_allchansn.mat'),'data_avg_allchans');
% save(fullfile(config{1}.datasavedir,'All_rats_data_avg_chantoplot.mat'),'data_avg_chantoplot');
% save(fullfile(config{1}.datasavedir,'All_rats_data_TFR.mat'),'data_TFR');
% save(fullfile(config{1}.datasavedir,'All_rats_Seizure_Infos.mat'),'Seizure_Infos');


%% Figures all rats
for do_normalize = [true false]
    
    [config] = dtx_setparams_probe([]);
    
    % Same datasavedir for all rats
    load(fullfile(config{1}.datasavedir,'All_rats_data_avg_allchansn.mat'),'data_avg_allchans');
    load(fullfile(config{1}.datasavedir,'All_rats_data_avg_chantoplot.mat'),'data_avg_chantoplot');
    load(fullfile(config{1}.datasavedir,'All_rats_data_TFR.mat'),'data_TFR');
    load(fullfile(config{1}.datasavedir,'All_rats_Seizure_Infos.mat'),'Seizure_Infos');
    
    % keep same data structure as for rat by rat analysis
    ipart = 1;
    irat = 1; %get the config of the first rat and modify specific fields
    
    %set config
    config{irat}.prefix = 'AllRats-';
    config{irat}.imagesavedir = fullfile(config{irat}.imagesavedir,'..','AllPatients');
    
    
    % Append data : one rat is one trial
    for imarker = 1:length(config{irat}.LFP.name)
        
        isdata = [];
        isdata = ~cellfun('isempty',data_avg_allchans{imarker}); %remove empty arrays for append
        if any(isdata)
            data_avg_allchans_allrats{ipart}{imarker}                = ft_appenddata([],data_avg_allchans{imarker}{isdata});
        else
            data_avg_allchans_allrats{ipart}{imarker}                = [];
        end
        
        isdata = [];
        isdata = ~cellfun('isempty',data_avg_chantoplot{imarker}); %remove empty arrays for append
        if any(isdata)
            data_avg_chantoplot_allrats{ipart}{imarker}                  = ft_appenddata([],data_avg_chantoplot{imarker}{isdata});
        else
            data_avg_chantoplot_allrats{ipart}{imarker}                = [];
        end
        
    end
    
    
    %normalize by dividing by value of max SlowWave at t=0 (alignment).
    %normalize emg by dividing by max(EMG) beetween 0 and 2
    if do_normalize
        config{irat}.prefix = 'Allrats-normalized-';
        config{irat}.imagesavedir = [config{irat}.imagesavedir,'-normalized'];
        
        for imarker = 1:length(config{irat}.LFP.name)
            
            for itrial = 1:length(config) %one trial is one rat
                %set norm values
                %norm_eeg : value at t=0 for align channel for eeg, or for
                %c4c3 channel for emg
                %norm_emg : max of trial between 0 and 2s
                t = [];
                norm_eeg = 0;
                if ~isempty(data_avg_chantoplot_allrats{ipart}{imarker})
                    if itrial <= length(data_avg_chantoplot_allrats{ipart}{imarker}.trial)
                        t_debut = data_avg_chantoplot_allrats{ipart}{imarker}.time{itrial}>-0.5;
                        t_fin   = data_avg_chantoplot_allrats{ipart}{imarker}.time{itrial}<0.5;
                        t_norm  = logical(t_debut .* t_fin);
                        norm_eeg = max(data_avg_chantoplot_allrats{ipart}{imarker}.trial{itrial}(t_norm));
                    end
                end
                
                %apply norm values to data
                if ~isempty(data_avg_allchans_allrats{ipart}{imarker})
                    if itrial <= length(data_avg_allchans_allrats{ipart}{imarker}.trial)
                        data_avg_allchans_allrats{ipart}{imarker}.trial{itrial} = data_avg_allchans_allrats{ipart}{imarker}.trial{itrial} / norm_eeg;
                    end
                end
                
                if ~isempty(data_avg_chantoplot_allrats{ipart}{imarker})
                    if itrial <= length(data_avg_chantoplot_allrats{ipart}{imarker}.trial)
                        data_avg_chantoplot_allrats{ipart}{imarker}.trial{itrial} = data_avg_chantoplot_allrats{ipart}{imarker}.trial{itrial} / norm_eeg;
                    end
                end
                
            end
        end %imarker
        
    end %do_normalize
    
    
    %plot data : exactly the same protocol as rat by rat
    dataEEG = []; %for figure with right and left eeg data, without emg data
    iEEG = 0;
    
    for imarker = 1:length(config{irat}.LFP.name)
        
        electrodetoplot = config{irat}.LFP.name{imarker}; %electrode of interest has been renamed over rats by get_patient_data.m
        dtx_plot_timecourse_eeg_emg(config{irat}, data_avg_chantoplot_allrats, ipart, imarker, 'eeg',electrodetoplot,[-inf, inf], true);
        dtx_plot_timecourse_eeg_emg(config{irat}, data_avg_chantoplot_allrats, ipart, imarker, 'eeg',electrodetoplot,[-5, 25], true);
        dtx_plot_timecourse_eeg_emg(config{irat}, data_avg_chantoplot_allrats, ipart, imarker, 'eeg',electrodetoplot,[-5, 5], true);
        dtx_plot_TFR(config{irat}, data_TFR, ipart, 1, true, true);
    end
    %figure common for 2 imarkers :
    imarker = 1;
    
    data_avg_allchanstoplot_allrats{ipart}{imarker} = ft_appenddata([], data_avg_chantoplot_allrats{ipart}{:}); %for comparing LFP M1 and EEG M1
    dtx_plot_comparison_severaleeg(config{irat},data_avg_allchanstoplot_allrats,ipart,imarker,{config{irat}.LFP.name{:}}, [-2 2], true, true)
    
    dtx_plot_overdraw_allchannels(config{irat},data_avg_allchans_allrats, ipart, imarker, true);
    dtx_plot_avg_allchannels(config{irat},data_avg_allchans_allrats,ipart,imarker,config{irat}.epoch.toi{imarker},true);
    dtx_plot_avg_allchannels(config{irat},data_avg_allchans_allrats,ipart,imarker,[-2 2],false);
    
end %do_normalize
%AJOUTER CV CV2 FANO FACTOR DE CHAQUE rat SUR UN PLOT


%
% %% Spike analysis
% config = dtx_setparams_probe([]);
% cfg = config{1}
%
% irat = 1;
%
% [MuseStruct]                     = readMuseMarkers(config{irat}, true);
%
% % align Muse markers according to peaks and detect whether they contain artefacts
% [MuseStruct]                     = alignMuseMarkers(config{irat},MuseStruct, false);
%
% [MuseStruct]                     = dtx_remove_wrong_seizure(config{irat}, MuseStruct,true, true);
%
% %remove seizure from 5s after SlowWave to Crise_End
% [MuseStruct_without_seizures]    = addMuseBAD(MuseStruct, 'all', 'all', 'SlowWave', 'Crise_End', 'all', 2, 1);
%
%
% [sampleinfo] = writeSpykingCircus(config{irat}, MuseStruct_without_seizures, true, true);
% writeSpykingCircusParameters(config{irat})
%
%
% % read spike-clustering results, and epoch around events
% [SpikeRaw] = readSpykingCircus_SpikeRaw(config{irat},true,'all');
% %[SpikeRaw, SpikeTrials] = readSpykingCircus(config{irat}, MuseStruct, false, 1);
%
% % compute event-related changes of spike rates, and other stats
% [stats_smooth, stats_binned] = spikeratestatsEvents(config{irat}, SpikeRaw, SpikeTrials, true);

%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     for macro_id = [config{irat}.LFP.electrodeToPlot(1), config{irat}.LFP.electrodeToPlot(2)]
%         %% time frequency analysis
%         for i_t_ftimwin = [9, 20, 40]
%             foi_max=50;
%
%             cfgtemp                         = [];
%             cfgtemp.channel                 = macro_id;%'all'; %ichannel;
%             cfgtemp.method                  = 'mtmconvol'; %relative change. 'mtmfft' : frequency content
%             cfgtemp.output                  = 'pow';
%             cfgtemp.taper                   = 'hanning';
%             cfgtemp.pad                     = 'nextpow2';
%             cfgtemp.keeptrials              = 'yes'; %return average
%             cfgtemp.foi                     = 1:0.5:50;
%             cfgtemp.t_ftimwin               = i_t_ftimwin./cfgtemp.foi;
%             %cfgtemp.t_ftimwin               = 20./cfgtemp.foi;
%             %cfgtemp.t_ftimwin               = 40./cfgtemp.foi;
%             %cfgtemp.t_ftimwin               = ones(size(cfgtemp.foi))*0.5;
%
%             cfgtemp.toi                     = [-5:0.01:25];
%             TFR_macro                    = ft_freqanalysis(cfgtemp,dat_macro{1}{1});
%
%             TFR_macro_log = TFR_macro;
%             TFR_macro_log.powspctrm = log(TFR_macro.powspctrm);
%
%             % without baseline relchange
%                         fig = figure;
%                         subplot(2,1,1);
%                         ft_singleplotTFR([], TFR_macro);
%
%                         title(sprintf('%s%s : Frequency power over time',config{irat}.prefix,dat_macro{1}{1}.label{macro_id}));
%                         xlim([-5, 25]);
%                         xlabel('Time from Slow Wave (s)');
%                         ylabel('Frequency (Hz)');
%
%                         subplot(2,1,2);
%                         ft_singleplotTFR([], TFR_macro_log);
%                         title(sprintf('%s%s : Frequency log power over time',config{irat}.prefix,dat_macro{1}{1}.label{macro_id}));
%                         xlim([-5, 25]);
%                         xlabel('Time from Slow Wave (s)');
%                         ylabel('Frequency (Hz)');
%
%

%
%             %ft_singleplotTFR([], TFR_macro_log);
%             %for foi_max = [25, 50]
%             foi_max=50;
%             %Plot
%             fig_title = sprintf('%s%s : Relative change from Baseline (-4s to -2s)',config{irat}.prefix,dat_macro{1}{1}.label{macro_id});
%             fig = figure;
%             subplot(2,1,1);
%             cfgtemp               = [];
%             %                 cfgtemp.channel         = macro_id;
%             cfgtemp.baseline        = [-4, -2];
%             cfgtemp.baselinetype    = 'relchange';
%             %cfgtemp.colorbar        = 'no';
%             cfgtemp.colorbar        = 'yes';
%             cfgtemp.zlim            = 'maxabs';
%             cfgtemp.xlim            = [-5 25];
%             cfgtemp.ylim            = [0 foi_max];
%             cfgtemp.title           = fig_title;
%             cfgtemp.parameter       = 'powspctrm';
%             cfgtemp.colormap        = parula(5000);
%             cfgtemp.renderer        = 'painters';
%             ft_singleplotTFR(cfgtemp, TFR_macro);
%
%             %                 if adapt_colormap==1
%             %colormap_axis = caxis;
%             %caxis([0,colormap_axis(2)]);
%             %                 end
%             %              colormap_axis = caxis;
%             %              caxis([colormap_axis(1)/1.8,colormap_axis(2)*1.5]);
%             %              xlabel('Time from Slow Wave (s)');
%             %              ylabel('Frequency (Hz)');
%             %              c=colorbar;
%             %              c.Label.String = 'Power relative change';
%             %
%             subplot(2,1,2);
%             fig2_title = sprintf('%s%s : Relative change (log power) from Baseline (-4s to -2s)',config{irat}.prefix,dat_macro{1}{1}.label{macro_id});
%             cfgtemp               = [];
%             %cfgtemp.channel         = macro_id;
%             cfgtemp.baseline        = [-4, -2];
%             cfgtemp.baselinetype    = 'relchange';
%             %cfgtemp.colorbar        = 'no';
%             cfgtemp.colorbar        = 'yes';
%             cfgtemp.zlim            = 'maxabs';
%             cfgtemp.xlim            = [-5 25];
%             cfgtemp.ylim            = [0 foi_max];
%             cfgtemp.title           = fig2_title;
%             cfgtemp.parameter       = 'powspctrm';
%             cfgtemp.colormap        = parula(5000);
%             cfgtemp.renderer        = 'painters';
%             ft_singleplotTFR(cfgtemp,TFR_macro_log);
%             %
%             %                 if adapt_colormap==1
%             %                     colormap_axis = caxis;
%             %caxis([-log10(colormap_axis(2)),log10(colormap_axis(2))]);
%             %                 end
%             % xlim([-5, 5]);
%             %              xlabel('Time from Slow Wave (s)');
%             %              ylabel('Frequency (Hz)');
%             %              c=colorbar;
%             %              c.Label.String = 'Log of power relative change';
%
%             % print to file
%             fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
%             set(fig,'PaperOrientation','landscape');
%             set(fig,'PaperUnits','normalized');
%             set(fig,'PaperPosition', [0 0 1 1]);
%             print(fig, '-dpdf', fullfile(config{irat}.imagesavedir,'tfr_seizures',[config{irat}.prefix,dat_macro{1}{1}.label{macro_id}, '_Avg_TFR_seizures_',num2str(foi_max),'Hz_',num2str(i_t_ftimwin)]),'-r600');
%             print(fig, '-dpng', fullfile(config{irat}.imagesavedir,'tfr_seizures',[config{irat}.prefix,dat_macro{1}{1}.label{macro_id}, '_Avg_TFR_seizures_',num2str(foi_max),'Hz_',num2str(i_t_ftimwin)]),'-r600');
%             close all
%
%             clear TFR_macro TFR_macro_log
%         end
%
%
%     end
%
% end
%





%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     %% Choose which electrode to plot : timecourse and TFR
%     %for macro_id = [config{irat}.LFP.electrodeToPlot(1), config{irat}.LFP.electrodeToPlot(2)]
%     for macro_id = config{irat}.LFP.electrodeToPlot(2)
%         %% plot all trials of one macro
%         for abscisse_max = [2, 5, 10, 25]
%             for i_size_amplitude = 1:4
%                 for i_amplification = [2, 5, 10]
%                     %for i_size_amplitude = 5:7
%                     fig             = figure;
%                     hold;
%                     fig.Renderer    = 'Painters'; % Else pdf is saved to bitmap
%                     h               = 25000/(10*(i_size_amplitude));
%                     n               = size(dat_macro{1}{1}.trial,2);
%                     for itrial = 1 : size(dat_macro{1}{1}.trial,2)
%                         plot(dat_macro{1}{1}.time{itrial},dat_macro{1}{1}.trial{itrial}(macro_id,:)*i_amplification + n*h - itrial*h,'k'); %first on top
%                     end
%
%                     xlabel('Time from SlowWave (s)');
%                     ylabel('Number of seizures');
%                     title(sprintf('%s : Aligned data from %s',config{irat}.prefix(1,1:end-1), dat_macro{1}{1}.label{macro_id}));
%                     set(gca, 'YTickLabel', '');
%                     tick = h;
%                     yticks(0 : tick*10 : n*h);
%                     yticklabels(n : -10 : 0);
%                     set(gca,'TickDir','out');
%                     axis tight
%                     if abscisse_max == 2
%                         xlim([-2, 2]);
%                     else
%                         xlim([-5, abscisse_max]);
%                     end
%
%                     % print to file
%                     set(fig,'PaperOrientation','landscape');
%                     set(fig,'PaperUnits','normalized');
%                     set(fig,'PaperPosition', [0 0 1 1]);
%                     set(fig,'Renderer','Painters');
%                     print(fig, '-dpdf', fullfile(config{irat}.imagesavedir,'timecourse_seizures',[config{irat}.prefix,dat_macro{1}{1}.label{macro_id},'_Timecourse_seizures_',...
%                         num2str(abscisse_max),'s_Sizex',num2str(i_amplification),'_',num2str(i_size_amplitude)]),'-r600');
%                     print(fig, '-dpng', fullfile(config{irat}.imagesavedir,'timecourse_seizures',[config{irat}.prefix,dat_macro{1}{1}.label{macro_id},'_Timecourse_seizures_',...
%                         num2str(abscisse_max),'s_Sizex',num2str(i_amplification),'_',num2str(i_size_amplitude)]),'-r600');
%                     close all
%                 end
%             end %i_size_amplitude
%         end %abscisse_max
%     end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %% average over trials for plotting of all channels
%     for abscisse_max = [2, 5, 10]
%
%         cfgtemp                 = [];
%         cfgtemp.vartrllength    = 2;
%         %dat_micro_rptavg        = ft_timelockanalysis(cfgtemp,dat_micro{1}{1});
%         dat_macro_rptavg        = ft_timelockanalysis(cfgtemp,dat_macro{1}{1});
%
%         fig=figure;
%         hold;
%         h=3000;
%         for i = 1 : size(dat_macro_rptavg.label,1)
%             %subplot(size(dat_macro_rptavg.label,1),1,size(dat_macro_rptavg.label,1)-i+1);
%             %         plot(dat_micro_rptavg.time,dat_micro_rptavg.avg)
%             plot(dat_macro_rptavg.time,dat_macro_rptavg.avg(i,:)+h*i,'color','black');
%
%         end
%         title(sprintf('Average of all SlowWaves of %s',config{irat}.prefix(1:end-1)));
%         axis tight;
%         xlabel('Time (s)');
%         ylabel('Channel name');
%         tick = h;
%         yticks(0 : tick : length(dat_macro_rptavg.label)*h);
%         set(gca, 'YTickLabel',[dat_macro_rptavg.label(end)', dat_macro_rptavg.label(1:end-1)']); %why Y axis not in the good order
%         set(gca,'TickDir','out');
%         if abscisse_max == 2
%             xlim([-2, 2]);
%         else
%             xlim([-5, abscisse_max]);
%         end
%         line ([abscisse_max-0.5,abscisse_max-0.5],[h*i-2500,h*i-500],'color','r','LineWidth',2);
%         text(abscisse_max-0.45, h*i-1400,'2mV','FontSize',15);
%
%
%         % print to file
%         fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
%         set(fig,'PaperOrientation','landscape');
%         set(fig,'PaperUnits','normalized');
%         set(fig,'PaperPosition', [0 0 1 1]);
%         print(fig, '-dpdf', fullfile(config{irat}.imagesavedir,'average_SlowWave',[config{irat}.prefix,'_Avg_SlowWaves_MACRO_',num2str(abscisse_max),'s']),'-r600');
%         print(fig, '-dpng', fullfile(config{irat}.imagesavedir,'average_SlowWave',[config{irat}.prefix,'_Avg_SlowWaves_MACRO_',num2str(abscisse_max),'s']),'-r600');
%
%         close all
%
%
%     end
%
% end
