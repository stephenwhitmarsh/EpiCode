function dtx_eegrodents_tfr(ipatient, config_script)

%config_script : name of the configuration script which will be called with
%'eval'. So this script can be used with several parameters scripts

%% Set parameters
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

%remove fieldtrip's output
ft_warning off
ft_notice off
ft_info off
ft_debug off
global ft_default
ft_default.checkconfig = 'silent';
ft_default.checkpath = 'once';
ft_default.showcallinfo = 'no';
ft_default.trackcallinfo = 'no';
ft_default.tracktimeinfo = 'no';

feature('DefaultCharacterSet', 'CP1252'); % To fix bug for weird character problems in reading neurlynx

config = eval(config_script);%dtx_setparams_eegvideo;

if ipatient > 0
    
    %% load precomputed data
    fname = fullfile(config{ipatient}.datasavedir,sprintf('%sLFP_cleanedforTFR.mat',config{ipatient}.prefix));
    if exist(fname,'file')
        fprintf('Reading %s\n',fname);
        temp = load(fname, 'data_cleaned');
        LFP{1}  = temp.data_cleaned;
    else
        fprintf('Reading %s\n', fullfile(config{ipatient}.datasavedir,[config{ipatient}.prefix,'LFP.mat']));
        temp = load(fullfile(config{ipatient}.datasavedir,[config{ipatient}.prefix,'LFP.mat']), 'LFP');
        LFP  = temp.LFP;
    end
    clear temp
    
    %% do time frequency analysis
    %setup baseline
    LFP{1}.Baseline = LFP{1}.SlowWave_begin;
    config{ipatient}.TFR.toi.Baseline = [-7 -2];
    
    %go trough each marker
    for markername = ["Baseline", "SlowWave_begin", "SlowWave", "Crise_End"]%string(config{ipatient}.TFR.name)
        
        fprintf('For marker %s\n', markername);
        
        %select only motor cortex electrode
        cfgtemp         = [];
        cfgtemp.channel = config{ipatient}.LFP.motorcortex;
        data_sel        = ft_selectdata(cfgtemp,LFP{1}.(markername));
        
        data_sel_avg    = ft_timelockanalysis([],data_sel);
        
        %high pass filter to avoid weird effect due to offset
        fprintf('patient %d : high pass filtering before TFR analysis\n', ipatient);
        cfgtemp             = [];
        cfgtemp.hpfilter    = 'yes';
        cfgtemp.hpfreq      = config{ipatient}.TFR.foi(1);
        cfgtemp.hpfilttype  = 'fir';
        data_sel            = ft_preprocessing(cfgtemp,data_sel);
        
        %exception for some rats because different way of putting crise_end
        %markers
        if strcmp(markername, "Crise_End")
            if strfind(config{ipatient}.prefix, 'Rat') %eeg video awake
                for itrial = 1:size(data_sel.trial,2)
                    data_sel.time{itrial} = data_sel.time{itrial} - 0.3;
                end
            end
            if ismember(config{ipatient}.prefix, {'DTX2-', 'DTX4-', 'DTX5-','DTX6-', 'DTX7-'})
                for itrial = 1:size(data_sel.trial,2)
                    data_sel.time{itrial} = data_sel.time{itrial} - 0.5;
                end
            end
            if ismember(config{ipatient}.prefix, {'DTX52-'})
                for itrial = 1:size(data_sel.trial,2)
                    data_sel.time{itrial} = data_sel.time{itrial} - 1;
                end
            end
        end
       
        %compute TFR
        fprintf('patient %d : do TFR analysis\n', ipatient);
             
        cfgtemp            = [];
        cfgtemp.channel    = 'all';
        cfgtemp.method     = 'mtmconvol';
        cfgtemp.output     = 'pow';
        cfgtemp.taper      = 'dpss';
        cfgtemp.tapsmofrq  = config{ipatient}.TFR.tapsmofrq; %number, the amount of spectral smoothing through multi-tapering. Note that 4 Hz smoothing means plus-minus 4 Hz, i.e. a 8 Hz smoothing box.
        cfgtemp.pad        = 'nextpow2';
        cfgtemp.keeptrials = 'no';
        cfgtemp.foi        = config{ipatient}.TFR.foi;
        cfgtemp.t_ftimwin  = ones(size(cfgtemp.foi)) .* config{ipatient}.TFR.timewinsize; %setparam
        cfgtemp.toi        = config{ipatient}.TFR.toi.(markername)(1) - 2 : config{ipatient}.TFR.timewinstep : config{ipatient}.TFR.toi.(markername)(2) + 2;
        TFR.(markername).raw   = ft_freqanalysis(cfgtemp,data_sel);
        clear data_sel
        
        %select latency (remove pad)
        cfgtemp         = [];
        cfgtemp.latency = config{ipatient}.TFR.toi.(markername);
        TFR.(markername).raw = ft_selectdata(cfgtemp,TFR.(markername).raw);
        
        %remove cfg field as it uses lot of memory
        TFR.(markername).raw = rmfield(TFR.(markername).raw,'cfg');

        %check the consistency of dimord
        if ~strcmp(TFR.(markername).raw.dimord, 'chan_freq_time')%~strcmp(TFR.(markername).raw.dimord, 'rpt_chan_freq_time')
            error('powspctrm dimensions are not consistent between conditions');
        end
        
        %compute and correct baseline
        marker_baseline  = "Baseline";
        for idata = "raw"%["raw", "log"]
            cfgtemp          = [];
            cfgtemp.latency  = config{ipatient}.TFR.toi.Baseline;
            TFR.baseline.(idata) = ft_selectdata(cfgtemp, TFR.(marker_baseline).(idata));
            
            TFR.(markername).(sprintf('%s_blcorrected',idata)) = TFR.(markername).(idata);
            for ifreq = 1:size(TFR.(markername).raw.freq,2)
                baselinetemp = squeeze(nanmean(TFR.baseline.(idata).powspctrm(1,ifreq,:)));
                TFR.(markername).(sprintf('%s_blcorrected',idata)).powspctrm(1,ifreq,:) = (TFR.(markername).(idata).powspctrm(1,ifreq,:) - baselinetemp) ./ baselinetemp;
            end
        end
        
        %plot data
        for i_blcorrected = ["","_blcorrected"]
            fig = figure;
            sgtitle(sprintf('%s %s : %s (%d trials)',config{ipatient}.prefix(1:end-1),TFR.(markername).raw.label{1}, markername,size(LFP{1}.(markername).trial,2)),'Interpreter','none','Fontsize',20);
            for idata = "raw"%["raw","log"]
                                 
%                 baseline = TFR.Baseline.(sprintf('%s%s',idata,i_blcorrected));
%                  
%                 %reduce baseline duration for smaller markers
%                 duration = TFR.(markername).(sprintf('%s%s',idata,i_blcorrected)).time(end) - TFR.(markername).(sprintf('%s%s',idata,i_blcorrected)).time(1);
%                 cfgtemp        = [];
%                 cfgtemp.latency = [baseline.time(end)-duration, baseline.time(end)];
%                 baseline   = ft_selectdata(cfgtemp,baseline);
%                 
%                 if strcmp(i_blcorrected,"_blcorrected")
%                     baseline_zeros = baseline;
%                     baseline_zeros.powspctrm = zeros(size(baseline.powspctrm));
%                 end
%                                 
%                 %% compute cluster-based permutation statistics
%                 
%                 cfgtemp = [];
%                 cfgtemp.channel          = 1;
%                 cfgtemp.latency          = 'all';
%                 cfgtemp.frequency        = 'all';
%                 cfgtemp.method           = 'montecarlo';
%                 cfgtemp.statistic        = 'ft_statfun_depsamplesT';
%                 cfgtemp.correctm         = 'cluster';
%                 cfgtemp.clusteralpha     = 0.05;
%                 cfgtemp.clusterstatistic = 'maxsum';
%                 cfgtemp.minnbchan        = 1;
%                 cfgtemp.tail             = 0;
%                 cfgtemp.clustertail      = 0;
%                 cfgtemp.alpha            = 0.025;
%                 cfgtemp.numrandomization = 1000;
%                 nb_rpt = size(baseline.powspctrm,1);
%                 cfgtemp.design(1,1:nb_rpt)          = 1:nb_rpt;%rpt
%                 cfgtemp.design(1,nb_rpt+1:2*nb_rpt) = 1:nb_rpt;%rpt
%                 cfgtemp.design(2,1:nb_rpt)          = 1;
%                 cfgtemp.design(2,nb_rpt+1:2*nb_rpt) = 2;
%                 cfgtemp.uvar     = 1;
%                 cfgtemp.ivar     = 2;
%                 
%                 stats.(sprintf('%s%s',idata,i_blcorrected)).baseline_all   = ft_freqstatistics(cfgtemp, TFR.(markername).(sprintf('%s%s',idata,i_blcorrected)), baseline);
%                 stats.(sprintf('%s%s',idata,i_blcorrected)).baseline_all   = rmfield(stats.(sprintf('%s%s',idata,i_blcorrected)).baseline_all, 'cfg');
%                 stats.(sprintf('%s%s',idata,i_blcorrected)).baseline_zeros = ft_freqstatistics(cfgtemp, TFR.(markername).(sprintf('%s%s',idata,i_blcorrected)), baseline);
%                 stats.(sprintf('%s%s',idata,i_blcorrected)).baseline_zeros = rmfield(stats.(sprintf('%s%s',idata,i_blcorrected)).baseline_zeros, 'cfg');
%                 clear baseline baseline_all
                               
                
%                 %iplot = iplot+1;
%                 %subplot(2,1,iplot);hold on;
%                 TFR.(markername).(sprintf('%s%s',idata,i_blcorrected)).powspctrm = nanmean(TFR.(markername).(sprintf('%s%s',idata,i_blcorrected)).powspctrm,1);
%                 TFR.(markername).(sprintf('%s%s',idata,i_blcorrected)).powspctrm = permute(TFR.(markername).(sprintf('%s%s',idata,i_blcorrected)).powspctrm,[2 3 4 1]);
%                 TFR.(markername).(sprintf('%s%s',idata,i_blcorrected)).dimord = 'chan_freq_time';
%                 
                %plot TFR
                cfgtemp         = [];
                cfgtemp.zlim    = 'maxmin';
                cfgtemp.figure  = 'gcf';
                ft_singleplotTFR(cfgtemp, TFR.(markername).(sprintf('%s%s',idata,i_blcorrected)));
                
                ylabel('Frequency (Hz)');
                set(gca,'TickDir','out','FontWeight','bold', 'FontSize', 15);
                xlabel('Time (s)');
                title(sprintf('%s%s',idata,i_blcorrected),'Interpreter','none');
                xlim(config{ipatient}.TFR.toi.(markername));
                
                %plot avg LFP
                if contains(markername,'SlowWave')
                    yyaxis right
                    plot(data_sel_avg.time,data_sel_avg.avg,'white','LineWidth',2);
                    ylabel('uV');
                    yyaxis left
                end

                if strcmp(i_blcorrected, '_blcorrected')
                    switch idata
                        case "raw"
                            caxis([-1 1]);
                        case "log"
                            caxis([-0.5 0.5]);
                    end
                else
                    colormap_axis = caxis;
                    caxis([0 colormap_axis(2)]);
                end
            end
            if strcmp(i_blcorrected, '_blcorrected')
                subdir = 'bl_corrected';
            else
                subdir = 'raw';
            end
            fname = fullfile(config{ipatient}.imagesavedir,'..','TFR',subdir,char(markername),sprintf('%sTFR_%s%s',config{ipatient}.prefix,markername,i_blcorrected));
            dtx_savefigure(fig,fname,'png','pdf','close');
        end
    end
    
    %save TFR data
    fname = fullfile(config{ipatient}.datasavedir,sprintf('%sTFR.mat',config{ipatient}.prefix));
    fprintf('saving data to %s\n',fullfile(config{ipatient}.datasavedir,sprintf('%sTFR.mat',config{ipatient}.prefix)));
    save(fname,'TFR','-v7.3');
%     fname = fullfile(config{ipatient}.datasavedir,sprintf('%sTFR_stats.mat',config{ipatient}.prefix));
%     fprintf('saving stats to %s\n',fname);
%     save(fname,'stats','-v7.3');
    
elseif ipatient == 0 %grand average over patients
    for markername = ["Baseline", "SlowWave_begin", "SlowWave", "Crise_End"]%string(config{ipatient}.TFR.name)
    i=0;
    %load precomputed TFR and LFP
    for ipatient = 1:size(config,2)
        %         if contains(markername, 'SlowWave')
        %             if ismember(ipatient,[2 4 5 7 13])
        %                 continue
        %             end
        %         elseif contains(markername, 'Crise')
        %             if ismember(ipatient,[1 4])
        %                 continue
        %             end
        %         end
        i=i+1;
        %read TFR
        fname           = fullfile(config{ipatient}.datasavedir,sprintf('%sTFR.mat',config{ipatient}.prefix));
        fprintf('Reading %s\n',fname);
        temp            = load(fname,'TFR');
        TFR{i}   = temp.TFR;
        clear temp
        %read LFP
        fname               = fullfile(config{ipatient}.datasavedir,sprintf('%sLFP_cleanedforTFR.mat',config{ipatient}.prefix));
        fprintf('Reading %s\n',fname);
        temp                = load(fname, 'data_cleaned');
        temp    = temp.data_cleaned;
        for markername = ["SlowWave", "SlowWave_begin"]
            LFP_avg.(markername){i} = ft_timelockanalysis([],temp.(markername));
            LFP_avg.(markername){i}.label = {'M1G'};
        end
        clear temp
    end
    for markername = ["SlowWave", "SlowWave_begin"]
        cfgtemp = [];
        cfgtemp.latency = config{ipatient}.TFR.toi.(markername);
        LFP_grandavg.(markername) = ft_timelockgrandaverage(cfgtemp,LFP_avg.(markername){:});
    end
    
    %go trough each marker
    %for markername = ["Baseline", "SlowWave_begin", "SlowWave", "Crise_End"]%string(config{ipatient}.TFR.name)
        for i_blcorrected = ["","_blcorrected"]
            
            fprintf('For marker %s\n', markername);
                
            for idata = "raw"%,"log"]
                
                %reorganize data and baseline, to use it with ft_freqgrangaverage
                clear data
                for ipatient = 1:size(TFR,2)
                    %data
                    TFR{ipatient}.(markername).(sprintf('%s%s',idata,i_blcorrected)).label = {'M1G'};
                    data{ipatient} = TFR{ipatient}.(markername).(sprintf('%s%s',idata,i_blcorrected));
                    %baseline
                    %                     TFR{ipatient}.Baseline.(sprintf('%s%s',idata,i_blcorrected)).label = {'M1G'};
                    %                     baseline{ipatient} = TFR{ipatient}.Baseline.(sprintf('%s%s',idata,i_blcorrected));
                end
                cfgtemp = [];
                cfgtemp.keepindividual = 'no';
                TFR_avg.(markername).(sprintf('%s%s',idata,i_blcorrected)) = ft_freqgrandaverage(cfgtemp,data{:});
                
                %reduce baseline duration for smaller markers
                %                 duration = TFR_avg.(markername).(sprintf('%s%s',idata,i_blcorrected)).time(end) - TFR_avg.(markername).(sprintf('%s%s',idata,i_blcorrected)).time(1);
                %                 cfgtemp        = [];
                %                 cfgtemp.toilim = [baseline{1}.time(end)-duration, baseline{1}.time(end)];
                %                 cfgtemp.keepindividual = 'yes';
                %                 baseline_all   = ft_freqgrandaverage(cfgtemp,baseline{:});
                %
                %                 if strcmp(i_blcorrected,"_blcorrected")
                %                     baseline_zeros = baseline_all;
                %                     baseline_zeros.powspctrm = zeros(size(baseline_all.powspctrm));
                %                 end
                %
                %                 %% compute cluster-based permutation statistics
                %
                %                 cfgtemp = [];
                %                 cfgtemp.channel          = {'M1G'};
                %                 cfgtemp.latency          = 'all';
                %                 cfgtemp.frequency        = 'all';
                %                 cfgtemp.method           = 'montecarlo';
                %                 cfgtemp.statistic        = 'ft_statfun_depsamplesT';
                %                 cfgtemp.correctm         = 'cluster';
                %                 cfgtemp.clusteralpha     = 0.05;
                %                 cfgtemp.clusterstatistic = 'maxsum';
                %                 cfgtemp.minnbchan        = 1;
                %                 cfgtemp.tail             = 0;
                %                 cfgtemp.clustertail      = 0;
                %                 cfgtemp.alpha            = 0.025;
                %                 cfgtemp.numrandomization = 1000;
                %                 cfgtemp.design(1,1:size(data,2))                = 1:size(data,2);
                %                 cfgtemp.design(1,size(data,2)+1:2*size(data,2)) = 1:size(data,2);
                %                 cfgtemp.design(2,1:size(data,2))                = 1;
                %                 cfgtemp.design(2,size(data,2)+1:2*size(data,2)) = 2;
                %                 cfgtemp.uvar     = 1;
                %                 cfgtemp.ivar     = 2;
                %
                %                 stats.baseline_all   = ft_freqstatistics(cfgtemp, TFR_avg.(markername).(sprintf('%s%s',idata,i_blcorrected)), baseline_all);
                %                 stats.baseline_zeros = ft_freqstatistics(cfgtemp, TFR_avg.(markername).(sprintf('%s%s',idata,i_blcorrected)), baseline_zeros);
                %                 clear baseline baseline_all
                %
                %
                %                 %% plot grand averaged TFR
                %                 cfgtemp = [];
                %                 cfgtemp.keepindividual = 'no';
                %                 TFR_avg.(markername).(sprintf('%s%s',idata,i_blcorrected)) = ft_freqgrandaverage([],data{:});
                %
                %iplot = iplot+1;
                %subplot(2,1,iplot);hold on;
                
                for i_maxfreq = [200,100,50]
                    
                    fig = figure;
                    sgtitle(sprintf('%s : grandaverage over patients',markername),'Interpreter','none','Fontsize',20);
                    
                    
                    cfgtemp = [];
                    cfgtemp.zlim = 'maxabs';
                    cfgtemp.figure = 'gcf';
                    cfgtemp.xlim = config{ipatient}.TFR.toi.(markername);
                    cfgtemp.ylim = [config{ipatient}.TFR.foi(1), i_maxfreq];
                    ft_singleplotTFR(cfgtemp, TFR_avg.(markername).(sprintf('%s%s',idata,i_blcorrected)));
                    ylabel('Frequency (Hz)');
                    xlim(config{ipatient}.TFR.toi.(markername));
                    ylim([config{ipatient}.TFR.foi(1), i_maxfreq]);
                    xlabel('Time (s)');
                    set(gca,'TickDir','out','FontWeight','bold', 'FontSize', 15);
                    
                    %plot avg LFP
                    if contains(markername,'SlowWave')
                        yyaxis right
                        plot(LFP_grandavg.(markername).time,LFP_grandavg.(markername).avg,'white','LineWidth',2);
                        ylabel('uV');
                        yyaxis left
                    end
                    
                    if strcmp(i_blcorrected, '_blcorrected')
                        %colormap_axis = caxis;
                        switch idata
                            case "raw"
                                caxis([-1 1]);
                            case "log"
                                caxis([-0.5 0.5]);
                        end
                    end
                    
                    if strcmp(i_blcorrected, '_blcorrected')
                        subdir = 'bl_corrected';
                    else
                        subdir = 'raw';
                    end
                    fname = fullfile(config{ipatient}.imagesavedir,'..','TFR',subdir,char(markername),sprintf('Allpatients-TFR_%s%s_%d',markername,i_blcorrected,i_maxfreq));
                    dtx_savefigure(fig,fname,'png','pdf','fig','close'); %save to fig to plot lfp average afterwards if needed
                end
            end
        end
    end
    %save TFR average and stats
    %     save(fullfile(config{ipatient}.datasavedir, 'allpatients_TFR_stats.mat'), 'stats', '-v7.3');
    save(fullfile(config{ipatient}.datasavedir, 'allpatients_TFR_avg.mat'), 'TFR_avg', '-v7.3');
end


end


