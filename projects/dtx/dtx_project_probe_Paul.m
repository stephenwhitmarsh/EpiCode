
%% Analyse of DTX seizures


% Stephen Whitmarsh, modified by Paul Baudin
%
% requires bandpassFilter.m from Mario
% requires releaseDec2015 from Neuralynx website
%requires signal processing toolbox


addpath C:\Users\paul.baudin\Documents\MATLAB\fieldtrip;
addpath (genpath('C:\Users\paul.baudin\Documents\MATLAB\EpiCode\projects\dtx'));
addpath (genpath('C:\Users\paul.baudin\Documents\MATLAB\EpiCode\external'));
addpath (genpath('C:\Users\paul.baudin\Documents\MATLAB\EpiCode\shared'));
addpath C:\Users\paul.baudin\Documents\MATLAB\DTX;
addpath C:\Users\paul.baudin\Documents\MATLAB\MatlabImportExport_v6.0.0;
ft_defaults



feature('DefaultCharacterSet', 'CP1252') % To fix bug for weird character problems in reading neurlynx





config = dtx_setparams_probe([]);
cfg = config{1}

%% LFP analysis
    
for irat = 1:6
    %% Get right LFP data
    % read muse markers
    [MuseStruct]    = readMuseMarkers(config{irat}, true);
    
    % align Muse markers according to peaks and detect whether they contain artefacts
    [MuseStruct]    = alignMuseMarkers(config{irat},MuseStruct, true);
    
    [MuseStruct] =   dtx_remove_wrong_seizure(config{irat}, MuseStruct, true);
    [dat_LFP] = readLFP(config{irat}, MuseStruct, false, false);
    %dat_LFP{ipart}{imarker}

    ipart = 1;
    imarker = 1;
    
    if irat == 2
        electrodeToPlot = 'ECoGS1';
    else
        electrodeToPlot = 'ECoGM1G';
    end
    
    dtx_plot_timecourse_eeg_emg(config{irat}, dat_LFP, ipart, imarker, 'eeg',electrodeToPlot,[-inf, inf], true);
    dtx_plot_timecourse_eeg_emg(config{irat}, dat_LFP, ipart, imarker, 'eeg',electrodeToPlot,[-5, 25], true);

    %AJOUTER DANS CONFIG SI ON VEUT QUE LE COMPTAGE DES CRISES PRENNE EN
    %COMPTE LE TEMPS ENTRE DEUX FICHIERS OU L'IGNORE (car �a peut �tre un
    %pb par exemple si seulement 1 ou 0 crise par fichier)
    
    %En fonction du r�sultat de l'alignement : voi si mrk EMG � placer
    %plut�t par rapport � OL ou par rapport � d�but EMG. 
    %Ce serait mieux par rappor t� d�but EMG : permet de moyenner en
    %fonction de la r�ponse EMG. Et possibilit� si besoin de choisir le LFP
    %correspondant pour au contraire montrer la r�ponse par rapport au LFP
    %Et pour les patients : Replacer les marqueurs EMG du coup
   
    
    
    
    
end
 dtx_plot_count_seizure(config{irat}, MuseStruct, ipart, true); %need to be done with all markers, not markers removed wrong seizures


%% Spike analysis
config = dtx_setparams_probe([]);
cfg = config{1}

irat = 1;

[MuseStruct]                     = readMuseMarkers(config{irat}, true);

% align Muse markers according to peaks and detect whether they contain artefacts
[MuseStruct]                     = alignMuseMarkers(config{irat},MuseStruct, false);

[MuseStruct]                     = dtx_remove_wrong_seizure(config{irat}, MuseStruct,true, true);

%remove seizure from 5s after SlowWave to Crise_End
[MuseStruct_without_seizures]    = addMuseBAD(MuseStruct, 'all', 'all', 'SlowWave', 'Crise_End', 'all', 2, 1);


[sampleinfo] = writeSpykingCircus(config{irat}, MuseStruct_without_seizures, true, true);
writeSpykingCircusParameters(config{irat})


% read spike-clustering results, and epoch around events
[SpikeRaw] = readSpykingCircus_SpikeRaw(config{irat},true,'all');
%[SpikeRaw, SpikeTrials] = readSpykingCircus(config{irat}, MuseStruct, false, 1);

% compute event-related changes of spike rates, and other stats
[stats_smooth, stats_binned] = spikeratestatsEvents(config{irat}, SpikeRaw, SpikeTrials, true);

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
