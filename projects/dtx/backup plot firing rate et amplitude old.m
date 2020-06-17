   %% Plot firing rate along all data
%             if hasevent, subplot(7,4,1:3);hold; else, subplot(4,4,1:3);hold; end
%             
%             %for baseline
%             for itrial = 1:size(SpikeTrials{ipart}{baseline_index}.trialinfo,1)
%                 yyaxis left
%                 clear spike_indx
%                 spike_indx              = SpikeTrials{ipart}{baseline_index}.trial{i_unit} == SpikeTrials{ipart}{baseline_index}.trialinfo(itrial,7);
%                 trial_begin             = SpikeTrials{ipart}{baseline_index}.trialinfo(itrial, 3) / spike_Fs;
%                 trial_end               = SpikeTrials{ipart}{baseline_index}.trialinfo(itrial, 4) / spike_Fs;
%                 trial_freq_avg(itrial)  = 1 / nanmean(stats{ipart}.(cfg.name{baseline_index}).isi.isi{i_unit}(spike_indx));
%                 baseline              = plot([trial_begin trial_end], log10([trial_freq_avg(itrial) trial_freq_avg(itrial)]),'-b','LineWidth',2);
%                 
%                 yyaxis right
%                 if log10(trial_freq_avg(itrial)) ~= Inf && log10(trial_freq_avg(itrial)) ~= -Inf && ~isnan(trial_freq_avg(itrial))
%                     hasdata = plot([trial_begin trial_end],[5 5],'-k','LineWidth', 2);
%                 end
%             end
%             stats{ipart}.(cfg.name{baseline_index}).freq_trialavg{i_unit} = trial_freq_avg;
%             
%             %for event
%             if hasevent
%                 clear trial_freq_avg
%                 for itrial = 1:size(SpikeTrials{ipart}{ievent}.trialinfo,1)
%                     yyaxis left
%                     clear spike_indx
%                     spike_indx              = SpikeTrials{ipart}{ievent}.trial{i_unit} == SpikeTrials{ipart}{ievent}.trialinfo(itrial,7);
%                     trial_begin             = SpikeTrials{ipart}{ievent}.trialinfo(itrial, 3) / spike_Fs;
%                     trial_end               = SpikeTrials{ipart}{ievent}.trialinfo(itrial, 4) / spike_Fs;
%                     trial_freq_avg(itrial)  = 1 / nanmean(stats{ipart}.(cfg.name{ievent}).isi.isi{i_unit}(spike_indx));
%                     event                   = scatter(trial_begin, log10(trial_freq_avg(itrial)),7,'or','filled');%, 'filled');
%                     yyaxis right
%                     if log10(trial_freq_avg(itrial)) ~= Inf && log10(trial_freq_avg(itrial)) ~= -Inf && ~isnan(trial_freq_avg(itrial))
%                         hasdata = plot([trial_begin trial_end],[5 5],'-k','LineWidth', 2);
%                     end
%                 end
%                 stats{ipart}.(cfg.name{ievent}).freq_trialavg{i_unit} = trial_freq_avg;
%             end
%             
%             yyaxis left
%             axis tight
%             ax = axis;
%             xticks(0:3600:ax(2));
%             xticklabels(xticks/3600); %convert seconds to hours
%             xlim([0 Inf]);
%             xlabel('Time (hours)');
%             set(gca, 'YColor', 'k');
%             ylabel(sprintf('Mean firing rate \nof each trial'));
%             yticklabels(10.^yticks);
%             setfig();
%             
%             yyaxis right
%             ylim([0 100]);
%             set(gca, 'YColor', 'none');
%             
%             if exist('hasdata', 'var')
%                 if hasevent, legend([baseline, event, hasdata],cfg.name{baseline_index},cfg.name{ievent},'Has data','location','eastoutside'); end
%                 if ~hasevent, legend([baseline, hasdata],cfg.name{baseline_index},'Has data','location','eastoutside'); end
%             end
%             
            %% Plot amplitude along all data %FIXME : do it after solving the issue of selecting amplitudes in trials
%             if hasevent, subplot(7,4,5:7);hold; else, subplot(4,4,5:7);hold; end
%             
%             clear baseline event hasdata
%             
%             %for baseline
%             for itrial = 1:size(SpikeTrials{ipart}{baseline_index}.trialinfo,1)
%                 clear spike_indx
%                 spike_indx              = SpikeTrials{ipart}{baseline_index}.trial{i_unit} == SpikeTrials{ipart}{baseline_index}.trialinfo(itrial,7);
%                 trial_begin             = SpikeTrials{ipart}{baseline_index}.trialinfo(itrial, 3) / spike_Fs;
%                 trial_end               = SpikeTrials{ipart}{baseline_index}.trialinfo(itrial, 4) / spike_Fs;
%                 trial_amp_avg(itrial)   = nanmean(SpikeTrials{ipart}{baseline_index}.amplitude{i_unit}(spike_indx));
%                 trial_amp_std(itrial)   = nanstd(SpikeTrials{ipart}{baseline_index}.amplitude{i_unit}(spike_indx));
%                 baseline                = plot([trial_begin trial_end], [trial_amp_avg(itrial) trial_amp_avg(itrial)],'-b','LineWidth',2);
%             end
%             stats{ipart}.(cfg.name{baseline_index}).amp_trialavg.avg{i_unit} = trial_amp_avg;
%             stats{ipart}.(cfg.name{baseline_index}).amp_trialavg.std{i_unit} = trial_amp_std;
%             
%             %for event
%             if hasevent
%                 clear trial_amp_avg trial_amp_avg
%                 for itrial = 1:size(SpikeTrials{ipart}{ievent}.trialinfo,1)
%                     clear spike_indx
%                     spike_indx              = SpikeTrials{ipart}{ievent}.trial{i_unit} == SpikeTrials{ipart}{ievent}.trialinfo(itrial,7);
%                     trial_begin             = SpikeTrials{ipart}{ievent}.trialinfo(itrial, 3) / spike_Fs;
%                     trial_end               = SpikeTrials{ipart}{ievent}.trialinfo(itrial, 4) / spike_Fs;
%                     trial_amp_avg(itrial)   = nanmean(SpikeTrials{ipart}{ievent}.amplitude{i_unit}(spike_indx));
%                     trial_amp_std(itrial)   = nanstd(SpikeTrials{ipart}{ievent}.amplitude{i_unit}(spike_indx));
%                     event                   = scatter(trial_begin, trial_amp_avg(itrial),7,'or','filled');
%                     
%                 end
%                 stats{ipart}.(cfg.name{ievent}).amp_trialavg.avg{i_unit} = trial_amp_avg;
%                 stats{ipart}.(cfg.name{ievent}).amp_trialavg.std{i_unit} = trial_amp_std;
%             end
%             
%             axis tight
%             ax = axis;
%             xticks(0:3600:ax(2));
%             xticklabels(xticks/3600); %convert seconds to hours
%             xlim([0 Inf]);
%             xlabel('Time (hours)');
%             set(gca, 'YColor', 'k');
%             ylabel(sprintf('Mean spike amplitude \nof each trial (µV)'));
%             setfig();
%             
%             if hasevent,  legend([baseline, event],cfg.name{baseline_index},cfg.name{ievent},'location','eastoutside'); end
%             if ~hasevent, legend(baseline,cfg.name{baseline_index},'location','eastoutside'); end
%             
%             yyaxis right %to set the same size of plot as firing rata
%             set(gca, 'YColor', 'none');