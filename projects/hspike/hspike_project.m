%% Analysis script for Spike x Sleep x SUA project
%
% (c) Stephen Whitmarsh, stephen.whitmarsh@gmail.com
%

%% Add path

if isunix
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/projects/hspike/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/shared/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/fieldtrip
    addpath(genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/scripts/releaseDec2015/'));
end

if ispc
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\projects\hspike
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\shared
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\fieldtrip
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\MatlabImportExport_v6.0.0 % to read neuralynx files faster
end

ft_defaults

feature('DefaultCharacterSet', 'CP1252') % To fix bug for weird character problems in reading neurlynx


%% TODO: 
% add saving of results in plotHypnogramStats

%% General analyses, looping over patients

for ipatient = 1
    
    % load settings
    config = hspike_setparams;
    
    % export hypnogram to muse
    export_hypnogram(config{ipatient});
    
    % read muse markers
    [MuseStruct] = readMuseMarkers(config{ipatient}, false);
    
    % align Muse markers according to peaks
    [MuseStruct] = alignMuseMarkers(config{ipatient},MuseStruct, false);
    
    % automatically detect spikes
    detectSpikes(config{ipatient}, MuseStruct, true, true);
    
    % read hypnogram as table
    [PSGtable] = PSG2table(config{ipatient}, MuseStruct, true);
    
    % plot hypnogram
    plotHypnogram(config{ipatient},MuseStruct);
    
    % events vs. hypnogram statistics and plots
    [MuseStruct, marker, hypnogram] = hypnogramStats(config{ipatient}, MuseStruct, true);
     
    % calculate TFR over all files, per part
    TFR = doTFRcontinuous(config{ipatient}, MuseStruct, true);
    
    % plot TFR 
    plotTFRcontinuous(config{ipatient},TFR);
    
    % write data concatinated for SC, artefacts, and output sampleinfo per file
    writeSpykingCircus(config{ipatient}, MuseStruct, false, false);
    
    % write parameters for spyking circus
    writeSpykingCircusParameters(config{ipatient})
    
    % read spike-clustering results, and epoch around events
    [SpikeRaw, SpikeTrials] = readSpykingCircus(config{ipatient}, MuseStruct, false, 'all');
    
    % compute event-related changes of spike rates, and other stats
    [stats_smooth, stats_binned] = spikeratestatsEvents(config{ipatient}, SpikeRaw, SpikeTrials, false);
    
    % read spike-clustering results, and label according to polysomnograpy
    [SpikeRawPSG, SpikeTrialsPSG] = readSpykingCircusPSG(config{ipatient}, MuseStruct, false, 'all');
    
    % computer spike stats, and label according to polysomnography
    [SpikeStatsPSG] = spikeratestatsPSG(config{ipatient}, SpikeRawPSG, SpikeTrialsPSG, hypnogram, true);
    
    % read LFP data
    [LFP] = readLFP(config{ipatient}, MuseStruct, true, true);
    
     
    TFRlog = TFR;
    TFRlog.powspctrm = TFRlog.powspctr
    % plot TFR
    

    
    
    
    
    
    
    % append nights
    cfg = [];
    cfg.keepsampleinfo = 'no';
    dat_micro_append{1} = ft_appenddata(cfg,dat_micro{1}{1},dat_micro{2}{1},dat_micro{3}{1});
    dat_macro_append{1} = ft_appenddata(cfg,dat_macro{1}{1},dat_macro{2}{1},dat_macro{3}{1});
    dat_micro_append{2} = ft_appenddata(cfg,dat_micro{1}{2},dat_micro{2}{2},dat_micro{3}{2});
    dat_macro_append{2} = ft_appenddata(cfg,dat_macro{1}{2},dat_macro{2}{2},dat_macro{3}{2});
    
    % average
    dat_micro_avg{1} = ft_timelockanalysis([],dat_micro_append{1});
    dat_macro_avg{1} = ft_timelockanalysis([],dat_macro_append{1});
    dat_micro_avg{2} = ft_timelockanalysis([],dat_micro_append{2});
    dat_macro_avg{2} = ft_timelockanalysis([],dat_macro_append{2});
    
    % plot LFP
    fig = figure;
    
    subplot(2,2,1); hold;
    plot(dat_macro_avg{1}.time, dat_macro_avg{1}.avg','k')
    ylim([-650,500]);
    y = ylim;
    plot([0 0],[y(1), y(2)],'k:');
    title('Manual Spike Detection Macro');
    axis tight
    xlabel('time (s)');
    ylabel('uVolts');
    
    subplot(2,2,3); hold;
    plot(dat_micro_avg{1}.time, dat_micro_avg{1}.avg','k')
    ylim([-650,500]);
    y = ylim;
    plot([0 0],[y(1), y(2)],'k:');
    title('Manual Spike Detection Macro');
    axis tight
    xlabel('time (s)');
    ylabel('uVolts');
    
    subplot(2,2,2); hold;
    plot(dat_macro_avg{2}.time, dat_macro_avg{2}.avg','k')
    ylim([-650,500]);
    y = ylim;
    plot([0 0],[y(1), y(2)],'k:');
    title('Automatic Spike Detection Macro');
    axis tight
    xlabel('time (s)');
    ylabel('uVolts');
    
    subplot(2,2,4); hold;
    plot(dat_micro_avg{2}.time, dat_micro_avg{2}.avg','k')
    ylim([-650,500]);
    y = ylim;
    plot([0 0],[y(1), y(2)],'k:');
    title('Automatic Spike Detection Macro');
    axis tight
    xlabel('time (s)');
    ylabel('uVolts');
    
    % print ISI to file
    fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
    set(fig,'PaperOrientation','landscape');
    set(fig,'PaperUnits','normalized');
    set(fig,'PaperPosition', [0 0 1 1]);
    print(fig, '-dpdf', fullfile(config{ipatient}.imagesavedir,'averageLFP.pdf'),'-r300');
    
    % TFR
    % time frequency analysis
    cfgtemp                         = [];
    cfgtemp.channel                 = 'all'; %ichannel;
    cfgtemp.method                  = 'mtmconvol';
    cfgtemp.output                  = 'pow';
    cfgtemp.taper                   = 'hanning';
    cfgtemp.pad                     = 'nextpow2';
    cfgtemp.keeptrials              = 'yes';
    cfgtemp.foi                     = 1:150;
    cfgtemp.t_ftimwin               = 7./cfgtemp.foi;
    cfgtemp.toi                     = config{ipatient}.epoch.toi{1}(1):0.01: config{ipatient}.epoch.toi{1}(2);
    TFR_micro{1}                    = ft_freqanalysis(cfgtemp,dat_micro_append{1});
    TFR_micro{2}                    = ft_freqanalysis(cfgtemp,dat_micro_append{2});
    TFR_macro{1}                    = ft_freqanalysis(cfgtemp,dat_macro_append{1});
    TFR_macro{2}                    = ft_freqanalysis(cfgtemp,dat_macro_append{2});
    
    save(fullfile(config{ipatient}.imagesavedir,[config{ipatient}.prefix, 'TFR']),'TFR*','-v7.3');
    load(fullfile(config{ipatient}.imagesavedir,[config{ipatient}.prefix, 'TFR']));
    
    % plot TFR
    
    fig = figure;
    cfgtemp               = [];
    cfgtemp.channel         = 'all';
    cfgtemp.baseline        = [-0.1, 0];
    %     cfgtemp.baselinetype    = 'relchange';
    cfgtemp.colorbar        = 'no';
    cfgtemp.colorbar        = 'yes';
    cfgtemp.zlim            = 'maxabs';
    %     cfgtemp.xlim            = config{ipatient}.;
    %     cfgtemp.title           = 'Relative change from Baseline';
    cfgtemp.parameter       = 'powspctrm';
    cfgtemp.colormap        = parula(5000);
    cfgtemp.renderer        = 'painters';
    ft_singleplotTFR(cfgtemp,TFR_micro{1});
    
    % print ISI to file
    fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
    set(fig,'PaperOrientation','landscape');
    set(fig,'PaperUnits','normalized');
    set(fig,'PaperPosition', [0 0 1 1]);
    print(fig, '-dpdf', fullfile(config{ipatient}.imagesavedir,[config{ipatient}.prefix, 'Avg_TFR_seizures.pdf']),'-r600');
    
    
    
end







%     dat_macro_append.trialinfo(dat_macro_append.trialinfo(:,4) == -1,4) = 0;

figure; hold;
for i = unique(dat_macro_append.trialinfo(:,4))'
    subplot(2,5,i+1);
    cfg = [];
    cfg.trials = find(dat_macro_append.trialinfo(:,4) == i);
    ft_singleplotER(cfg,dat_macro_append);
    
    %         subplot(2,5,i+1+6); hold;
    %         for ii = 1 : size(dat_macro_append.trialinfo,1)
    %             if dat_macro_append.trialinfo(ii,4) == i
    %                 plot(dat_macro_append.trial{ii});
    %             end
    %         end
end

% put all in none matrix
dat = zeros(size(dat_macro_append{1}.trial,2),size(dat_macro_append{1}.trial{1},2));
for i = 1 : size(dat_macro_append{1}.trial,2)
    dat(i,:) = dat_macro_append{1}.trial{i}(1,:);
    cls(i) = dat_macro_append{1}.trialinfo(i,4);
end

% combine for classifier
d = [dat, cls'];

% for neural net

cls_dummy = zeros(size(dat,1),size(unique(cls),2));
for i = 1 : size(dat,1)
    cls_dummy(i,cls(i)+1) = 1;
end

% cluster pattern
findPattern(config{ipatient}, dat_micro, dat_macro_append, 1)






% Example:
%   ini = IniConfig();
%   ini.ReadFile('example.ini')
%   sections = ini.GetSections()
%   [keys, count_keys] = ini.GetKeys(sections{1})
%   values = ini.GetValues(sections{1}, keys)
%   new_values(:) = {rand()};
%   ini.SetValues(sections{1}, keys, new_values, '%.3f')
%   ini.WriteFile('example1.ini')
%
% Example:
%   ini = IniConfig();
%   ini.AddSections({'Some Section 1', 'Some Section 2'})
%   ini.AddKeys('Some Section 1', {'some_key1', 'some_key2'}, {'hello!', [10, 20]})
%   ini.AddKeys('Some Section 2', 'some_key3', true)
%   ini.AddKeys('Some Section 2', 'some_key1')
%   ini.WriteFile('example2.ini')
%
% Example:
%   ini = IniConfig();
%   ini.AddSections('Some Section 1')
%   ini.AddKeys('Some Section 1', 'some_key1', 'hello!')
%   ini.AddKeys('Some Section 1', {'some_key2', 'some_key3'}, {[10, 20], [false, true]})
%   ini.WriteFile('example31.ini')
%   ini.RemoveKeys('Some Section 1', {'some_key1', 'some_key3'})
%   ini.RenameKeys('Some Section 1', 'some_key2', 'renamed_some_key2')
%   ini.RenameSections('Some Section 1', 'Renamed Section 1')
%   ini.WriteFile('example32.ini')
%






% average over trials for plotting
cfgtemp                 = [];
cfgtemp.vartrllength    = 2;
dat_micro_rptavg        = ft_timelockanalysis(cfgtemp,dat_micro{1});
dat_macro_rptavg        = ft_timelockanalysis(cfgtemp,dat_macro{1});


figure;
for i = 1 : size(dat_macro_rptavg.label,1)
    subplot(size(dat_macro_rptavg.label,1),1,i);
    %         plot(dat_micro_rptavg.time,dat_micro_rptavg.avg)
    plot(dat_macro_rptavg.time,dat_macro_rptavg.avg(i,:));
    title(dat_macro_rptavg.label{i});
end

figure;
for i = 1 : size(dat_micro_rptavg.label,1)
    subplot(size(dat_micro_rptavg.label,1),1,i);
    %         plot(dat_micro_rptavg.time,dat_micro_rptavg.avg)
    plot(dat_micro_rptavg.time,dat_micro_rptavg.avg(i,:));
    title(dat_micro_rptavg.label{i});
end


% select time interval
cfgtemp                 = [];
cfgtemp.latency         = [cfg.epoch.toi{imarker}(1) cfg.epoch.toi{imarker}(2)];
dat_micro_rptavg        = ft_selectdata(cfgtemp,dat_micro_rptavg);
dat_macro_rptavg        = ft_selectdata(cfgtemp,dat_macro_rptavg);







% plot LFP timecourse examples for article
% plotTimeCourses(config{ipatient});

% plot LFP data
[FFT_micro_trials,TFR_micro_trials,TFR_macro_trials,stat_TFR_micro] = plotLFP(config{ipatient}, dat_micro, dat_macro, true);

% write data concatinated for SC, and update config with sampleinfo
config{ipatient}                        = writeSpykingCircus(config{ipatient}, MuseStruct_micro, false);

% read raw spike data from SC, and segment into trials
[SpikeRaw, SpikeTrials]                 = readSpykingCircus(config{ipatient}, MuseStruct_micro, false);

% read and plot spikerate overview, and get the stats
[SpikeRateStats, stats_bar, sdf_orig_out, sdf_bar_out, corrs] = spikeratestats(config{ipatient}, SpikeRaw, SpikeTrials, false);

% read and plot LFP of spike events
% [spike_LFP, spike_LFP_avg]  = spikeLFP(config{ipatient},SpikeRaw);


%% plot correlations between micro and macro

for ipatient =  1 : 3
    
    config = setparams([]);
    
    % read muse markers
    [MuseStruct_micro, MuseStruct_macro]    = readMuseMarkers(config{ipatient}, false);
    
    % align Muse markers according to peaks and detect whether they contain artefacts
    [MuseStruct_micro, MuseStruct_macro]    = alignMuseMarkers(config{ipatient},MuseStruct_micro, MuseStruct_macro, false);
    
    % read LFP data
    [dat_micro, dat_macro] = readLFP(config{ipatient},MuseStruct_micro, MuseStruct_macro, false);
    
    % plot LFP timecourse examples for article
    % plotTimeCourses(config{ipatient});
    
    % plot LFP data
    [FFT_micro_trials{ipatient}, TFR_micro_trials{ipatient}, TFR_macro_trials{ipatient}, stat_TFR_micro{ipatient}, corrs{ipatient}] = plotLFP(config{ipatient}, dat_micro, dat_macro, false);
end

close all
fig = figure; hold;
mloc = [3, 2, 1; 3, 2, 1; 2, 1, 0];
xi = 0;
c = [];
labels = {};
for ipatient =  1 : 3
    for imarker = 1 : size(stat_TFR_micro{ipatient},2)
        xi = xi + 1
        
        for icontact = 1 : length(stat_TFR_micro{ipatient}{imarker}.corrs.avg)
            xpos = (xi-1) * 3
            ypos = (icontact-1) * 3
            c = stat_TFR_micro{ipatient}{mloc(ipatient,imarker)}.corrs.avg(icontact);
            if c >= 0
                col = 'g.';
            else
                col = 'r.';
            end
            plot(xpos,ypos,col,'markersize',abs(c*300));
            labels{xi} = config{ipatient}.name{mloc(ipatient,imarker)};
            text(xpos,ypos+1.5,sprintf('%.2f',c),'HorizontalAlignment','center');
        end
    end
end
xticks((0:7)*3);
xticklabels(labels);
yticks((0:6)*3);
yticklabels([1:7]);
xlim([-3 8*3]);
ylim([-3 7*3]);
xlabel('Pattern');
ylabel('Macro Contact');

% print ISI to file
fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
set(fig,'PaperOrientation','landscape');
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition', [0 0 1 1]);
print(fig, '-dpdf', '/network/lustre/iss01/charpier/stephen.whitmarsh/WhitmarshEpilepsy/correlations_LFP_macro.pdf','-r300');





%% correlate spikerate with LFP, and plot firing-rates for all units
for ipatient = 1 : 3
    
    config = setparams([]);
    
    % read muse markers
    [MuseStruct_micro, MuseStruct_macro]    = readMuseMarkers(config{ipatient}, false);
    
    % align Muse markers according to peaks and detect whether they contain artefacts
    [MuseStruct_micro, MuseStruct_macro]    = alignMuseMarkers(config{ipatient},MuseStruct_micro, MuseStruct_macro, false);
    
    % read LFP data
    [dat_micro, dat_macro] = readLFP(config{ipatient},MuseStruct_micro, MuseStruct_macro, false);
    
    % read raw spike data from SC, and segment into trials
    [SpikeRaw, SpikeTrials]                 = readSpykingCircus(config{ipatient}, MuseStruct_micro, false);
    
    % read and plot spikerate overview, and get the stats
    [SpikeRateStats, stats_bar, sdf_orig_out, sdf_bar_out] = spikeratestats(config{ipatient}, SpikeRaw, SpikeTrials, false);
    
    for imarker = 1 : size(dat_micro,2)
        
        % find datafilename corresponding to channel
        channelindx = [];
        for ilabel = 1 : size(dat_micro{imarker}.label,1)
            if strfind(dat_micro{imarker}.label{ilabel},config{ipatient}.align.channel{imarker})
                channelindx = ilabel;
            end
        end
        
        cfg = [];
        cfg.channel = channelindx;
        dat_micro_sel = ft_timelockanalysis(cfg,dat_micro{imarker});
        
        cfg = [];
        cfg.latency = [sdf_orig_out{imarker}{1}.time(1), sdf_orig_out{imarker}{1}.time(end)];
        dat_micro_sel = ft_selectdata(cfg,dat_micro_sel);
        
        for itemp = 1 : size(sdf_orig_out{imarker},2)
            [corrs{ipatient}{imarker}{itemp}.rho, corrs{ipatient}{imarker}{itemp}.p] = corr(sdf_orig_out{imarker}{itemp}.avg(itemp,:)',dat_micro_sel.avg','rows','complete');
        end
    end
    
end

%% compare baseline firingrates and create table (SEPARATE)
clear stat_bl SpikeRateStats SpikeRateStats_bar SpikeRaw SpikeTrials tbl unit
tbl = table;
ntemp = 1;


for ipatient = 1 : 3
    
    config = setparams([]);
    
    % read muse markers
    [MuseStruct_micro, MuseStruct_macro]    = readMuseMarkers(config{ipatient}, false);
    
    % align Muse markers according to peaks and detect whether they contain artefacts
    [MuseStruct_micro, MuseStruct_macro]    = alignMuseMarkers(config{ipatient},MuseStruct_micro, MuseStruct_macro, false);
    
    config{ipatient}                        = writeSpykingCircus(config{ipatient}, MuseStruct_micro, false);
    
    % read raw spike data from SC, and segment into trials
    [SpikeRaw, SpikeTrials]                 = readSpykingCircus(config{ipatient}, MuseStruct_micro, false);
    
    % read and plot spikerate overview, and get the stats
    [SpikeRateStats{ipatient}, SpikeRateStats_bar{ipatient}, ~, ~]              = spikeratestats(config{ipatient}, SpikeRaw, SpikeTrials, false);
    
    % compare baselines
    for itemp = 1 : size(SpikeRateStats_bar{ipatient}.clusterstat{1},2)
        p = [];
        x = [];
        for ilabel = 1 : size(SpikeRateStats_bar{ipatient}.clusterstat,2)
            p = [p; ones(length(SpikeRateStats_bar{ipatient}.clusterstat{ilabel}{itemp}.bl.trialavg),1) * ilabel];
            x = [x; SpikeRateStats_bar{ipatient}.clusterstat{ilabel}{itemp}.bl.trialavg];
        end
        mdl = fitlm(x,p);
        stat_bl{ipatient}{itemp} = anova(mdl,'summary');
    end
    
    % create table
    for itemp = 1 : size(SpikeRateStats_bar{ipatient}.clusterstat{1},2)
        tbl.nodule(ntemp)   = ipatient;
        total               = length(SpikeRateStats{ipatient}.isi{itemp});
        short               = sum(SpikeRateStats{ipatient}.isi{itemp} <= 2);
        long                = sum(SpikeRateStats{ipatient}.isi{itemp} > 2);
        tbl.unit{ntemp}     = itemp;
        tbl.total(ntemp)    = total;
        
        if short/total * 100 < 1
            tbl.percRPV{ntemp} = sprintf('$%.2f^s$',round(short/total * 100,2));
        else
            tbl.percRPV{ntemp} = sprintf('$%.2f^m$',round(short/total * 100,2));
        end
        
        unit(ntemp)         = itemp;
        
        for ilabel = 1 : size(SpikeRateStats_bar{ipatient}.clusterstat,2)
            
            if isfield(SpikeRateStats_bar{ipatient}.clusterstat{ilabel}{itemp},'posclusters')
                if ~isempty(SpikeRateStats_bar{ipatient}.clusterstat{ilabel}{itemp}.posclusters)
                    if SpikeRateStats_bar{ipatient}.clusterstat{ilabel}{itemp}.posclusters(1).prob < 0.025 % stats are sorted to most significant
                        if SpikeRateStats_bar{ipatient}.clusterstat{ilabel}{itemp}.posclusters(1).prob < 0.001
                            tbl.([config{ipatient}.name{ilabel},'_increase']){ntemp} = sprintf('$%.0f^{xxx}$',SpikeRateStats_bar{ipatient}.clusterstat{ilabel}{itemp}.maxcluster.perc{1});
                        elseif SpikeRateStats_bar{ipatient}.clusterstat{ilabel}{itemp}.posclusters(1).prob < 0.01
                            tbl.([config{ipatient}.name{ilabel},'_increase']){ntemp} = sprintf('$%.0f^{xx}$', SpikeRateStats_bar{ipatient}.clusterstat{ilabel}{itemp}.maxcluster.perc{1});
                        else
                            tbl.([config{ipatient}.name{ilabel},'_increase']){ntemp} = sprintf('$%.0f^{x}$',  SpikeRateStats_bar{ipatient}.clusterstat{ilabel}{itemp}.maxcluster.perc{1});
                        end
                    else
                        tbl.([config{ipatient}.name{ilabel},'_increase']){ntemp} = 'n.s.';
                    end
                    
                else
                    tbl.([config{ipatient}.name{ilabel},'_increase']){ntemp} = 'n.s.';
                end
            else
                tbl.([config{ipatient}.name{ilabel},'_increase']){ntemp} = 'n.s.';
            end
            
            if isfield(SpikeRateStats_bar{ipatient}.clusterstat{ilabel}{itemp},'negclusters')
                if ~isempty(SpikeRateStats_bar{ipatient}.clusterstat{ilabel}{itemp}.negclusters)
                    if SpikeRateStats_bar{ipatient}.clusterstat{ilabel}{itemp}.negclusters(1).prob < 0.025 % stats are sorted to most significant
                        if SpikeRateStats_bar{ipatient}.clusterstat{ilabel}{itemp}.negclusters(1).prob < 0.001
                            tbl.([config{ipatient}.name{ilabel},'_decrease']){ntemp} = sprintf('$%.0f^{xxx}$',SpikeRateStats_bar{ipatient}.clusterstat{ilabel}{itemp}.mincluster.perc{1});
                        elseif SpikeRateStats_bar{ipatient}.clusterstat{ilabel}{itemp}.negclusters(1).prob < 0.01
                            tbl.([config{ipatient}.name{ilabel},'_decrease']){ntemp} = sprintf('$%.0f^{xx}$',SpikeRateStats_bar{ipatient}.clusterstat{ilabel}{itemp}.mincluster.perc{1});
                        else
                            tbl.([config{ipatient}.name{ilabel},'_decrease']){ntemp} = sprintf('$%.0f^{x}$',SpikeRateStats_bar{ipatient}.clusterstat{ilabel}{itemp}.mincluster.perc{1});
                        end
                    else
                        tbl.([config{ipatient}.name{ilabel},'_decrease']){ntemp} = 'n.s.';
                    end
                else
                    tbl.([config{ipatient}.name{ilabel},'_decrease']){ntemp} = 'n.s.';
                end
            else
                tbl.([config{ipatient}.name{ilabel},'_decrease']){ntemp} = 'n.s.';
            end
            
            
            if corrs{ipatient}{ilabel}{itemp}.p < 0.001
                tbl.([config{ipatient}.name{ilabel},'_corr']){ntemp}  = sprintf('$%.2f^{xxx}$',corrs{ipatient}{ilabel}{itemp}.rho);
            elseif corrs{ipatient}{ilabel}{itemp}.p < 0.01
                tbl.([config{ipatient}.name{ilabel},'_corr']){ntemp}  = sprintf('$%.2f^{xx}$',corrs{ipatient}{ilabel}{itemp}.rho);
            elseif corrs{ipatient}{ilabel}{itemp}.p < 0.05
                tbl.([config{ipatient}.name{ilabel},'_corr']){ntemp}  = sprintf('$%.2f^{x}$',corrs{ipatient}{ilabel}{itemp}.rho);
            else
                tbl.([config{ipatient}.name{ilabel},'_corr']){ntemp}  = 'n.s';
            end
            
            %             tbl.([config{ipatient}.name{ilabel},'_corr']) =
            
        end
        
        % zero crossing
        zci                 = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);
        
        % peak width accoridng to Gast et. al
        temp                = dir(fullfile(config{ipatient}.datasavedir,[config{ipatient}.prefix,'all_data_',config{ipatient}.circus.channel{1}(1:end-2),'_*.ncs']));
        hdr_fname           = fullfile(temp(1).folder,temp(1).name);
        hdr                 = ft_read_header(hdr_fname); % take the first file to extract the header of the data
        tempsel             = squeeze(SpikeRaw.template(itemp,SpikeRaw.template_maxchan(itemp),:));
        temptime            = ((1:size(SpikeRaw.template,3))/hdr.Fs*1000)';
        
        % interpolate template
        temptime_int        = linspace(temptime(1),temptime(end),10000);
        tempsel_int         = pchip(temptime,tempsel,temptime_int);
        [Ypos,Xpos]         = findpeaks(tempsel_int,temptime_int,'NPeaks',1,'SortStr','descend');
        [Yneg,Xneg]         = findpeaks(-tempsel_int,temptime_int,'NPeaks',2,'SortStr','descend','MinPeakDistance',0.5);
        %         plot([Xpos,Xneg(1)],[Yneg(1)*0.3, Yneg(1)*0.3],'-o','Color',[1 0 0],'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]);
        %         plot([Xpos,Xneg(2)],-[Yneg(2)*0.3, Yneg(2)*0.3],'-o','Color',[0 0 1],'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]);
        
        %         tbl.tp{ntemp} = sprintf('%.0f',abs(Xpos-Xneg(1))*1000);
        tbl.pt{ntemp}       = sprintf('%.0f',abs(Xpos-Xneg(2))*1000);
        tp(ntemp)           = abs(Xpos-Xneg(1))*1000;
        pt(ntemp)           = abs(Xpos-Xneg(2))*1000;
        
        midline             = min(tempsel_int) + (max(tempsel_int) - min(tempsel_int)) / 2 ;
        indx                = zci(tempsel_int - midline);
        tbl.width{ntemp}    = sprintf('%.0f',diff(temptime_int(indx))*1000);
        w(ntemp)            = diff(temptime_int(indx))*1000;
        
        
        ntemp = ntemp + 1;
    end
end

tbl.displayed{tbl.nodule == 1 & unit' == 2} = 'A';
tbl.displayed{tbl.nodule == 2 & unit' == 8} = 'B';
tbl.displayed{tbl.nodule == 2 & unit' == 5} = 'C';
tbl.displayed{tbl.nodule == 3 & unit' == 7} = 'D';

% tbl = sortrows(tbl,{'nodule','SUA'},{'ascend','descend'});

writetable(tbl,'/network/lustre/iss01/charpier/stephen.whitmarsh/WhitmarshEpilepsy/unittables.xls')


%% Make figure like Bartho 2000
fig = figure; hold;

i = strcmp(tbl.displayed, 'A') | strcmp(tbl.displayed, 'B') | strcmp(tbl.displayed, 'C') | strcmp(tbl.displayed, 'D');
scatter(w(i)/1000,pt(i)/1000,120,'r','d','filled');
xlabel('Half-amplitude duration');
ylabel('Trough to peak time');

i = strcmp(tbl.SUA, 'MUA');
scatter(w(i)/1000,pt(i)/1000,60,'k','linewidth',2);

i = strcmp(tbl.SUA, 'SUA');
scatter(w(i)/1000,pt(i)/1000,60,'k','filled');

for i = 1 : size(w,2)
    text(w(i)/1000+0.003,pt(i)/1000+0.015,sprintf('%d-%d',tbl.nodule(i),tbl.unit(i)));
end

axis square

% print ISI to file
fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
set(fig,'PaperOrientation','landscape');
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition', [0 0 1 1]);
print(fig, '-dpdf', '/network/lustre/iss01/charpier/stephen.whitmarsh/WhitmarshEpilepsy/scatter_spikes_time.pdf','-r600');


%% find max FFT
force = false;
for ipatient = 1 : 3
    config = setparams([]);
    [~, TFR_micro_trials{ipatient},~,~] = plotLFP(config{ipatient}, [], [], false);
end


for ipatient = 1 : 3
    for imarker = 1 : size(FFT_micro_trials{ipatient},2)
        figure;
        
        plot(FFT_micro_trials{ipatient}{imarker}.powspctrm');
    end
end


%% PLOT ALL FIRINGRATES


% load data

for ipatient = 1 : 3
    
    config = setparams([]);
    
    % read muse markers
    [MuseStruct_micro, MuseStruct_macro]    = readMuseMarkers(config{ipatient}, false);
    
    % align Muse markers according to peaks and detect whether they contain artefacts
    [MuseStruct_micro, MuseStruct_macro]    = alignMuseMarkers(config{ipatient},MuseStruct_micro, MuseStruct_macro, false);
    
    config{ipatient}                        = writeSpykingCircus(config{ipatient}, MuseStruct_micro, false);
    
    % read raw spike data from SC, and segment into trials
    [SpikeRaw, SpikeTrials]                 = readSpykingCircus(config{ipatient}, MuseStruct_micro, false);
    
    % read and plot spikerate overview, and get the stats
    [~, stats_bar{ipatient}, sdf_orig{ipatient}, sdf_bar{ipatient}]        = spikeratestats(config{ipatient}, SpikeRaw, SpikeTrials, false);
end

% plotting

fig = figure;
orient(fig,'portrait');
mloc = [3, 2, 1; 3, 2, 1; 2, 1, 0];
h = size(sdf_bar{1}{1},2) + size(sdf_bar{2}{1},2) + size(sdf_bar{3}{1},2) + 15; % add spacing between patients
w = 3;
hcum = cumsum([0 size(sdf_bar{1}{1},2)+1, size(sdf_bar{2}{1},2)+1]);
iplot = 1;
hi = 1;
cmap = parula;
for ipatient = 1 : 3
    for imarker = 1 : 3
        hi = 1;
        if mloc(ipatient,imarker) > 0
            
            for itemp = 1 : size(sdf_bar{ipatient}{imarker},2)
                
                wi = imarker;
                iplot = (hcum(ipatient) + hi-1)*w+wi;
                subplottight(h,w,iplot); hold;
                axis tight
                bar(sdf_bar{ipatient}{mloc(ipatient,imarker)}{itemp}.time,sdf_bar{ipatient}{mloc(ipatient,imarker)}{itemp}.avg,1,'facecolor',[127/255,127/255,127/255],'edgecolor',[127/255,127/255,127/255]);
                
                if isfield(stats_bar{ipatient}.clusterstat{mloc(ipatient,imarker)}{itemp},'posclusters')
                    for ipos = 1 : size(stats_bar{ipatient}.clusterstat{mloc(ipatient,imarker)}{itemp}.posclusters,2)
                        if stats_bar{ipatient}.clusterstat{mloc(ipatient,imarker)}{itemp}.posclusters(ipos).prob < config{ipatient}.stats.alpha
                            lag = size(sdf_bar{ipatient}{mloc(ipatient,imarker)}{itemp}.avg,1) - size(stats_bar{ipatient}.clusterstat{mloc(ipatient,imarker)}{itemp}.mask,2);
                            
                            sel = find(stats_bar{ipatient}.clusterstat{mloc(ipatient,imarker)}{itemp}.posclusterslabelmat == ipos);
                            bar(stats_bar{ipatient}.clusterstat{mloc(ipatient,imarker)}{itemp}.time(sel),sdf_bar{ipatient}{mloc(ipatient,imarker)}{itemp}.avg(sel+lag),1,'facecolor',[252/255,187/255,62/255],'edgecolor',[252/255,187/255,62/255]);
                        end
                    end
                end
                
                if isfield(stats_bar{ipatient}.clusterstat{mloc(ipatient,imarker)}{itemp},'negclusters')
                    for ipos = 1 : size(stats_bar{ipatient}.clusterstat{mloc(ipatient,imarker)}{itemp}.negclusters,2)
                        if stats_bar{ipatient}.clusterstat{mloc(ipatient,imarker)}{itemp}.negclusters(ipos).prob < config{ipatient}.stats.alpha
                            lag = size(sdf_bar{ipatient}{mloc(ipatient,imarker)}{itemp}.avg,1) - size(stats_bar{ipatient}.clusterstat{mloc(ipatient,imarker)}{itemp}.mask,2);
                            
                            sel = find(stats_bar{ipatient}.clusterstat{mloc(ipatient,imarker)}{itemp}.negclusterslabelmat == ipos);
                            bar(stats_bar{ipatient}.clusterstat{mloc(ipatient,imarker)}{itemp}.time(sel),sdf_bar{ipatient}{mloc(ipatient,imarker)}{itemp}.avg(sel+lag),1,'facecolor',[70/255,93/255,250/255],'edgecolor',[70/255,93/255,250/255]);
                        end
                    end
                end
                
                xt = xticks;
                set(gca,'TickDir','out');
                if itemp ~= size(sdf_bar{ipatient}{imarker},2)
                    set(gca,'xtick',[]);
                end
                hi = hi + 1;
            end
        end
    end
end

% print ISI to file
fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
set(fig,'PaperOrientation','portrait');
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition', [0 0 1 1]);
print(fig, '-dpdf', '/network/lustre/iss01/charpier/stephen.whitmarsh/WhitmarshEpilepsy/firing_rates_all.pdf','-r600');


%% plot ISIs and template morphologies

% load data
clear SpikeRaw SpikeTrials
for ipatient = 1 : 3
    
    config = setparams([]);
    
    % read muse markers
    [MuseStruct_micro, MuseStruct_macro]    = readMuseMarkers(config{ipatient}, false);
    
    % align Muse markers according to peaks and detect whether they contain artefacts
    [MuseStruct_micro, MuseStruct_macro]    = alignMuseMarkers(config{ipatient},MuseStruct_micro, MuseStruct_macro, false);
    
    config{ipatient}                        = writeSpykingCircus(config{ipatient}, MuseStruct_micro, false);
    
    % read raw spike data from SC, and segment into trials
    [SpikeRaw{ipatient}, SpikeTrials{ipatient}]                 = readSpykingCircus(config{ipatient}, MuseStruct_micro, false);
    
end

% plot

fig = figure;
orient(fig,'portrait');
mloc = [3, 2, 1; 3, 2, 1; 2, 1, 0];
h = size(sdf_bar{1}{1},2) + size(sdf_bar{2}{1},2) + size(sdf_bar{3}{1},2) + 15; % add spacing between patients
w = 1;

hcum = cumsum([0 size(sdf_bar{1}{1},2)+1, size(sdf_bar{2}{1},2)+1]);
iplot = 1;
cmap = parula;

for ipatient = 1 : 3
    
    for itemp = 1 : size(sdf_bar{ipatient}{1},2)
        
        subplottight(h,2,iplot); hold;
        
        temp        = dir(fullfile(config{ipatient}.datasavedir,[config{ipatient}.prefix,'all_data_',config{ipatient}.circus.channel{1}(1:end-2),'_*.ncs']));
        hdr_fname   = fullfile(temp(1).folder,temp(1).name);
        hdr         = ft_read_header(hdr_fname); % take the first file to extract the header of the data
        tempsel     = squeeze(SpikeRaw{ipatient}.template(itemp,SpikeRaw{ipatient}.template_maxchan(itemp),:));
        temptime    = ((1:size(SpikeRaw{ipatient}.template,3))/hdr.Fs*1000)';
        
        % interpolate template
        temptime_int = linspace(temptime(1),temptime(end),10000);
        tempsel_int = pchip(temptime,tempsel,temptime_int);
        plot(temptime_int,tempsel_int,'k');
        
        % zero crossing
        zci                 = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);
        axis tight
        [Ypos,Xpos] = findpeaks( tempsel_int,temptime_int,'NPeaks',1,'SortStr','descend');
        [Yneg,Xneg] = findpeaks(-tempsel_int,temptime_int,'NPeaks',2,'SortStr','descend','MinPeakDistance',0.5);
        
        plot([Xpos,Xneg(1)], [Yneg(1)*0.3, Yneg(1)*0.3],'-o','Color',[0 0 0],'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0],'MarkerSize',4);
        plot([Xpos,Xneg(2)],-[Yneg(2)*0.3, Yneg(2)*0.3],'-o','Color',[0 0 0],'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0],'MarkerSize',4);
        
        %       x = (Xpos + Xneg(1))/2;
        %       y = Yneg(1)*0.1;
        %       text(x,y,sprintf('%.0fms',abs(Xpos-Xneg(1))*1000),'HorizontalAlignment','center');
        %       x = (Xpos + Xneg(2))/2;
        %       y = -Yneg(2)*0.1;
        %       text(x,y,sprintf('%.0fms',abs(Xpos-Xneg(2))*1000),'HorizontalAlignment','center');
        midline = min(tempsel_int) + (max(tempsel_int) - min(tempsel_int)) / 2 ;
        indx = zci(tempsel_int - midline);
        plot(temptime_int(indx),[midline, midline],'-o','Color',[0 0 0],'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0],'MarkerSize',4);
        x = sum(temptime_int(indx))/length(indx);
        y = midline*1.1;
        %       text(x,y,sprintf('%.0fms',(temptime_int(indx(2))-temptime_int(indx(1)))*1000),'HorizontalAlignment','center');
        iplot = iplot + 1;
        
        xt = xticks;
        set(gca,'TickDir','out');
        set(gca,'ytick',[]);
        
        if itemp ~= size(sdf_bar{ipatient}{1},2)
            set(gca,'xtick',[]);
        else
            xlabel('time (ms)');
        end
        
        % ISI
        subplottight(h,2,iplot); hold;
        isi = diff(SpikeRaw{ipatient}.samples{itemp}) / hdr.Fs * 1000;
        histogram(isi,'BinWidth',0.5,'BinLimits',[0,25],'FaceColor',[0,0,0],'EdgeColor',[0,0,0],'FaceAlpha',1);
        %       bar(stats.isi_1s.time*1000,stats.isi_1s.avg(itemp,:),1);
        %       xticks(stats.isi_1s.time*1000);
        %       xtickangle(90);
        axis tight
        
        xt = xticks;
        set(gca,'TickDir','out');
        if itemp ~= size(sdf_bar{ipatient}{1},2)
            set(gca,'xtick',[]);
        else
            xlabel('time (ms)');
        end
        
        iplot = iplot + 1;
        
    end
    
    iplot = iplot + 2;
    
end

% print ISI to file
fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
set(fig,'PaperOrientation','portrait');
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition', [0 0 1 1]);
print(fig, '-dpdf', '/network/lustre/iss01/charpier/stephen.whitmarsh/data/images/templates_all.pdf','-r600');


























% for idir = 1 : size(quality.dat_mad_cnt,2)
%     for imad = 1 : 10
%         for ifile = 1 : size(quality.dat_mad_cnt{idir},2)
%             d(i)        = q_norm(imad,idir,ifile);
%             mad(i)      = imad;
%             dir(i)      = idir;
%             channel(i)  = ifile;
%             i = i + 1;
%         end
%     end
% end
%
% data = [mad', dir', channel', d'];
% labels = {'MAD','Directory','Channel','Count'};
%
% figure;
% [h, ax] = plotmatrix(data)
%
% for i = 1 : 4
%     xlabel(ax(4,i),labels{i});
%     ylabel(ax(i,1),labels{i});
% end








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOW RUN SPYKING-CIRCUS ON *-all_data_*.ncs %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read and plot all data as one big file from SC
%     cfg                     = [];
%     cfg.prefix              = [num2str(ipatient),'-'];
%     cfg.force               = false;
%     cfg.startend            = startend{ipatient};
%     cfg.prestim             = prestim{ipatient}; % +pad{ipatient};
%     cfg.poststim            = poststim{ipatient}; % +pad{ipatient};
%     cfg.imagesavedir        = imagesavedir{ipatient};
%     cfg.datasavedir         = datasavedir{ipatient};
%     cfg.label               = label{ipatient};
%     cfg.channel             = channel{ipatient};
%     cfg.resamplefs          = resamplefs{ipatient};
%     cfg.timwin              = timwin(ipatient,:);
%     cfg.suffix              = '-1';
%     cfg.sampleinfo          = sampleinfo;
%     [SpikeRaw, SpikeTrials, stat_x] = readSpykingCircus_allmarkers_alldata(cfg,MuseStruct_micro);

%
%
% %% Read and plot ERF around spike events
% cfg = [];
% cfg.prefix                  = [num2str(ipatient),'-'];
% cfg.ipatient                = ipatient;
% cfg.fnames_ncs              = fnames_ncs;
% cfg.datasavedir             = datasavedir{ipatient};
% cfg.channel                 = channel{ipatient};
% cfg.imagesavedir            = imagesavedir{ipatient};
% [spike_LFP, spike_LFP_avg]  = spikeLFP(cfg,SpikeRaw);
%
% for imarker = 1 : length(label{ipatient})
%
%     %% read LFP data
%     cfg                     = [];
%     cfg.force               = false; % whether to force recalculations
%     cfg.startend            = startend{ipatient}(imarker,:);
%     cfg.label               = label{ipatient}{imarker};
%     cfg.prestim             = prestim{ipatient}(imarker)  + pad{ipatient};
%     cfg.poststim            = poststim{ipatient}(imarker) + pad{ipatient};
%     cfg.hpfilter            = hpfilter{ipatient}{imarker};
%     cfg.hpfreq              = hpfreq{ipatient}{imarker};
%     cfg.baseline            = 'no';
%     cfg.micro_labels        = micro_labels{ipatient,imarker};
%     cfg.macro_labels        = macro_labels{ipatient,imarker};
%     cfg.resamplefs          = resamplefs{ipatient};
%     cfg.datasavedir         = datasavedir{ipatient};
%     [dat_micro, dat_macro]  = readEEG(cfg,MuseStruct_micro,MuseStruct_macro); % rename to readLFP
%
%     %% remove artefacts based on RMS - make function
%     cfg = [];
%     cfg.prestim             = prestim{ipatient}(imarker); % +pad{ipatient};
%     cfg.channel             = channel{ipatient}{imarker};
%     %         [dat_micro,dat_macro,c] = removeartefacts_corr(cfg,dat_micro, dat_macro); % remove c, an d place in plotLFP function
%     [~,~,c] = removeartefacts_corr(cfg,dat_micro, dat_macro); % remove c, an d place in plotLFP function
%
%     %% plot EEG data
%     cfg                     = [];
%     cfg.force               = false;
%     cfg.label               = label{ipatient}{imarker};
%     cfg.prestim             = prestim{ipatient}(imarker);
%     cfg.poststim            = poststim{ipatient}(imarker);
%     cfg.slidestep           = slidestep{ipatient}(imarker);
%     [Y, I]                  = sort(c,'descend');
%     cfg.representtrials     = I(1:20); % do within function
%     cfg.resamplefs          = resamplefs{ipatient};
%     cfg.channel             = channel{ipatient}{imarker};
%     cfg.binsize             = 0.1;
%     cfg.datasavedir         = datasavedir{ipatient};
%     cfg.imagesavedir        = imagesavedir{ipatient};
%     plotdata_corr(cfg,dat_micro,dat_macro);
%
%     %% read data at original high samplerate
%     cfg                     = [];
%     cfg.force               = false;
%     cfg.startend            = startend{ipatient}(imarker,:);
%     cfg.prestim             = prestim{ipatient}(imarker); % +pad{ipatient};
%     cfg.poststim            = poststim{ipatient}(imarker); % +pad{ipatient};
%     cfg.datasavedir         = datasavedir{ipatient};
%     cfg.label               = label{ipatient}{imarker};
%     cfg.channel             = micro_labels{ipatient,imarker};
%     cfg.hpfilter            = 'no';
%     cfg.hpfreq              = 100;
%     dat_microFs             = readMicroFs(cfg,MuseStruct_micro);
%
%     %% analyse spiketriggered
%     cfg = [];
%     cfg.prestim             = prestim{ipatient}(imarker); % +pad{ipatient};
%     cfg.poststim            = poststim{ipatient}(imarker); % +pad{ipatient};
%     cfg.datasavedir         = datasavedir{ipatient};
%     cfg.imagesavedir        = imagesavedir{ipatient};
%     cfg.prefix              = [num2str(ipatient),'-'];
%     cfg.fnames_ncs          = fnames_ncs;
%     spiketriggeredplots(cfg,dat_microFs,SpikeTrials,SpikeRaw);
%
%     %         % correlate between channels over trials
%     %         for itrial = 1 : size(dat_microFs.trial,2)
%     %             disp(num2str(itrial));
%     %             c(itrial,:,:) = corr(dat_microFs.trial{itrial}');
%     %         end
%     %
%     %         c_avg = squeeze(nanmean(c,1));
%     %         c_std = squeeze(nanstd(c,1));
%     %         c_3std = c_std * 3;
%     %
%     %         fig = figure;
%     %         im = image(c_avg*255);
%     %         colormap(jet(255));
%     %         im.CDataMapping = 'scaled';
%     %         colormap jet
%     %         [H,P,CI,STATS] = ttest(c-0.5);
%     %
%     %         % print to file
%     %         set(fig,'PaperOrientation','landscape');
%     %         set(fig,'PaperUnits','normalized');
%     %         set(fig,'PaperPosition', [0 0 1 1]);
%     %         print(fig, '-dpdf', fullfile(imagesavedir{ipatient},'correlation.pdf'),'-r600');
%     %         set(fig,'PaperOrientation','portrait');
%     %         print(fig, '-dpng', fullfile(imagesavedir{ipatient},'correlation.png'),'-r600');
%     %
%     %         %% connectivity analysis
%     %         cfg = [];
%     %         cfg.method = 'mtmfft';
%     %         cfg.taper = 'hanning';
%     %         cfg.output = 'fourier';
%     %         cfg.foi = 1:100;
%     %         cfg.pad = 'nextpow2';
%     %         freq = ft_freqanalysis(cfg,dat_microFs);
%     %
%     %         cfg = [];
%     %         cfg.method = 'coh';
%     %         coh = ft_connectivityanalysis(cfg,freq);
%     %
%     %         fig = figure;
%     %         cfg = [];
%     %         cfg.parameter = 'cohspctrm';
%     %         cfg.zlim = [0 1];
%     %         cfg.xlim = [1 40];
%     %         ft_connectivityplot(cfg, coh);
%     %
%     %
%     %         % print to file
%     %         set(fig,'PaperOrientation','landscape');
%     %         set(fig,'PaperUnits','normalized');
%     %         set(fig,'PaperPosition', [0 0 1 1]);
%     %         print(fig, '-dpdf', fullfile(imagesavedir{ipatient},'coherence.pdf'),'-r600');
%     %         print(fig, '-dpng', fullfile(imagesavedir{ipatient},'coherence.png'),'-r600');
%
% end % imarker
%
%
%
%
%
%
% close all
%
%
%
% %
% %
% % %% test correlations between nodules in patient #2
% %
% % cfgtemp                     = [];
% % cfgtemp.force               = true;
% % cfgtemp.datasavedir         = datasavedir{ipatient};
% % cfgtemp.label               = label{ipatient}{imarker};
% % cfgtemp.startend            = startend{ipatient};
% % cfgtemp.channel             = ["_TNmi_1","_TNmi_2","_TNmi_3","_TNmi_4","_TNmi_5","_TNmi_6","_TNmi_7","_TNmi_8", "mCasd_1","mCasd_2","mCasd_3","mCasd_4","mCasd_5","mCasd_6","mCasd_7","mCasd_8"];
% % dat_microFs                 = readMicroFs(cfgtemp,MuseStruct_micro);
% %
% %
% %
% %
% %
% %
% %
% %
% %
% %
% %
% %
% %
% %
% %
% %
% %
% %
% % % ilabel = 1;
% % %
% % % cfgtemp                         = [];
% % % cfgtemp.trl                     = cfg.trialinfo_single{ilabel} + cfg.trialinfo_all(ilabel,1) - 1; % + cfg.prestim(ilabel) * hdr.Fs;
% % % cfgtemp.trl(:,3)                = -ones(size(cfgtemp.trl,1),1) * cfg.prestim(ilabel) * hdr.Fs;
% % % cfgtemp.trlunit                 = 'samples';
% % % cfgtemp.hdr                     = hdr;
% % % SpikeTrials                     = ft_spike_maketrials(cfgtemp,SpikeRaw);
% %
% %
% % cfgtemp                     = [];
% % cfgtemp.force               = false;
% % cfgtemp.datasavedir         = datasavedir{ipatient};
% % cfgtemp.label               = label{ipatient}{ilabel};
% % dat_microFs                 = readMicroFs(cfgtemp,MuseStruct_micro);
% %
% % prefix                      = [num2str(ipatient),'-'];
% %
% % % add header - better to do it in readMicroFs
% % temp                        = dir(fullfile(datasavedir{ipatient},[prefix,'all_concatinated_',channel{ipatient}{imarker}(1:end-2),'_*.ncs']));
% %
% % hdr                         = ft_read_header(fullfile(temp(1).folder,temp(1).name)); % take the first file to extract the header of the data
% %
% % dat_microFs.hdr = hdr;
% % data_all = ft_appendspike([],dat_microFs,SpikeTrials);
% %
% %
% % cfgtemp                             = [];
% % cfgtemp.timwin                      = [-0.25 0.25];
% % cfgtemp.spikechannel                = channel{ipatient}{ilabel}; % first unit
% % cfgtemp.channel                     = micro_labels{ipatient,ilabel}; % first four chans
% % cfgtemp.latency                     = [-prestim{ipatient}(ilabel) poststim{ipatient}(ilabel)];
% % staPost                             = ft_spiketriggeredaverage(cfg, SpikeTrials);
% %
% % dat_micro
% % %% write deadfile for Spyking Circus
% % % writeSpykingCircusDeadFile(MuseStruct_micro_aligned)
% % % writeSpykingCircusDeadFile_concatinated(MuseStruct_micro_aligned)
% %

