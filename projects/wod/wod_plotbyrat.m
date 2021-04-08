function wod_plotbyrat(rat_list, configscript)

try %en local
    scriptpath = matlab.desktop.editor.getActiveFilename;
catch %cluster
    scriptpath = mfilename('fullpath');
end

epicodepath = [fileparts(fileparts(fileparts(fileparts(scriptpath)))), filesep];

addpath (genpath([epicodepath,'development']))
addpath (genpath([epicodepath,'shared']))
addpath (genpath([epicodepath,'external']))
addpath (genpath([epicodepath,'templates']))
addpath (genpath([epicodepath,'projects', filesep, 'wod']))
addpath (genpath([epicodepath,'projects', filesep, 'dtx']))

if ispc
    addpath \\lexport\iss01.charpier\analyses\wod\fieldtrip-20200607
elseif isunix
    addpath /network/lustre/iss01/charpier/analyses/wod/fieldtrip-20200607
end

ft_defaults

config = eval(configscript);



ipart = 1; %ipart is always 1 for this project



analysis_names = {'timefreq_wod', 'timefreq_wod_timenorm', 'timefreq_recovery','timefreq_wor','timefreq_wor_timenorm', 'timefreq_baseline','timefreq_wod_blcorrected', 'timefreq_wod_timenorm_blcorrected', 'timefreq_recovery_blcorrected','timefreq_wor_blcorrected','timefreq_wor_timenorm_blcorrected', 'timefreq_baseline_blcorrected','log_timefreq_wod', 'log_timefreq_wod_timenorm', 'log_timefreq_recovery','log_timefreq_wor','log_timefreq_wor_timenorm', 'log_timefreq_baseline','log_timefreq_wod_blcorrected', 'log_timefreq_wod_timenorm_blcorrected', 'log_timefreq_recovery_blcorrected','log_timefreq_wor_blcorrected','log_timefreq_wor_timenorm_blcorrected', 'log_timefreq_baseline_blcorrected'};

%go through each data structure
for idata = 1:size( analysis_names,2)
%% TFR by chan by rat
    
    
    for irat = rat_list %size(config,2)
       
        %load data
        fprintf('Load timefreq data for rat %d/%d\n', irat,size(config,2));
        data_temp = load(fullfile(config{irat}.datasavedir,[config{irat}.prefix,analysis_names{idata},'.mat']));
        
        data_rat = [];
        count_trials                    = 0;
        %get data from one rat in a structure
        for itrial = 1:size(data_temp.(analysis_names{idata}),2)
            count_trials = count_trials +1;
            wod_rat(count_trials) = count_trials; %store rat id for each wod (some rats have several wod)
            chan_list                                   = fieldnames(data_temp.(analysis_names{idata}){itrial});

            data_rat= data_temp.(analysis_names{idata}){itrial};
            
            for ichan= 1: numel(chan_list)
                chan_name = chan_list{ichan};
                %replace empty channels by nans
                hasdata = find(~structfun(@isempty,data_rat));
                if isempty(hasdata)
                    continue
                end
                
                if isempty(data_rat.(chan_name))
                    data_rat.(chan_name)                       = data_rat.(chan_list{hasdata(1)});
                    data_rat.(chan_name).powspctrm             = ones(size(data_rat.(chan_name).powspctrm));
                end
                
                data_plot = data_rat.(chan_name);
                
               
                fig=figure;
                
                sgtitle([analysis_names{idata},chan_name,sprintf('Rat %d/%d',irat,size(config,2)),' ',sprintf('WOD %d/%d\n',wod_rat(count_trials),size(data_temp.(analysis_names{idata}),2))], 'Interpreter', 'none', 'FontWeight', 'bold', 'FontSize', 23);
                subplot(2,1,1)
                %plot TFR
                %voir les paramètres optionnels dans le descriptifs de la fonction pour
                %modifier l'aspect du TFR. Avec les paramètres par défaut :
                cfgtemp         = [];
                cfgtemp.channel = 'all';
                cfgtemp.interactive = 'no';
                cfgtemp.colormap= 'jet';
                cfg.fontsize = 12;
                cfgtemp.ylim= [51 100];
                cfgtemp.masknans    = 'yes';
                ft_singleplotTFR(cfgtemp, data_plot);
                ft_pimpplot(fig, jet(5000))
                
                
                title('Time-Frequency [50 : 100]');
                %figure settings (fontsize,fontweight, ticks dir, renderer etc.)
                set(gca, 'TickDir', 'out', 'FontSize', 12, 'FontWeight', 'bold');
                xlabel('Time (s)');
                ylabel(sprintf('Frequency (Hz)'));
                axis tight;
                
                subplot(2,1,2)
                %plot TFR
                %voir les paramètres optionnels dans le descriptifs de la fonction pour
                %modifier l'aspect du TFR. Avec les paramètres par défaut :
                cfgtemp         = [];
                cfgtemp.channel = 'all';
                cfgtemp.interactive = 'no';
                cfgtemp.masknans    = 'yes';
                cfgtemp.colormap= 'jet';
                cfgtemp.ylim= [0 50];
                cfg.fontsize = 12;
                ft_singleplotTFR(cfgtemp, data_plot);
                ft_pimpplot(fig, jet(5000))
                
                
                title('Time-Frequency [1 : 50]');
                %figure settings (fontsize,fontweight, ticks dir, renderer etc.)
                set(gca, 'TickDir', 'out', 'FontSize', 12, 'FontWeight', 'bold');
                xlabel('Time (s)');
                ylim([1 50]);
                ylabel(sprintf('Frequency (Hz)'));
                axis tight;
                
                if ~isfolder(config{4}.imagesavedir)
                    fprintf('Creating directory %s\n', config{4}.imagesavedir);
                    mkdir(config{4}.imagesavedir)
                end
                
                rat_imagedir= fullfile(config{4}.imagesavedir_data{1},'by_rat',sprintf(config{irat}.prefix));
                
                
                if idata >= 7
                    rat_imagedir = fullfile(rat_imagedir, 'baseline_corrected');
                end
                
                if idata > 12 & idata < 19
                    rat_imagedir = fullfile(rat_imagedir, 'Log');
                end
                
                if idata >= 19
                    rat_imagedir = fullfile(rat_imagedir, 'Log','baseline corrected');
                end
                                        
                
                
                if ~isfolder(rat_imagedir)
                    fprintf('Creating directory %s\n',  rat_imagedir);
                    mkdir(rat_imagedir)
                end
                
                
                
                set(fig,'PaperOrientation','landscape');
                set(fig,'PaperUnits','normalized');
                set(fig,'PaperPosition', [0 0 1 1]);
                print(fig, '-dpdf', fullfile(rat_imagedir,sprintf('Rat_%d_WoD_%d_%s_%s.pdf',irat,wod_rat(count_trials),analysis_names{idata}, chan_name)),'-r600');
                print(fig, '-dpng', fullfile(rat_imagedir,sprintf('Rat_%d_WoD_%d_%s_%s.png',irat,wod_rat(count_trials),analysis_names{idata}, chan_name)),'-r600');
                %savefig(fig,fullfile(rat_imagedir,sprintf('Rat_%d_WoD_%d_%s_%s.pdf',irat,wod_rat(count_trials),analysis_names{idata}, chan_name)));
                close all
                
                
                
                
                
                
            end %ichan
            
            
            
            
        end %itrial
        clear data_temp data_plot
    end %irat
end %idata
