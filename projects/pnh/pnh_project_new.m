%% Analysis script for Valerio Frazzini's study on electrophysiology of pathoanatomy
%
% (c) Stephen Whitmarsh, stephen.whitmarsh@gmail.com
%
% requires bandpassFilter.m from Mario
% requires releaseDec2015 from Neuralynx website

%% Add path

restoredefaultpath

if isunix
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/fieldtrip
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/projects/pnh/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/shared/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/shared/utilities
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/subaxis
    addpath(genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/scripts/releaseDec2015/binaries'));
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/altmany-export_fig-8b0ba13/
    ft_defaults
end

if ispc
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\fieldtrip
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\projects\pnh
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\shared
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\shared\utilities
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external\subaxis
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\MatlabImportExport_v6.0.0 
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external\altmany-export_fig-8b0ba13\
    ft_defaults
end

% initialise
MuseStruct{4}               = [];
SpikeRaw{4}                 = [];
SpikeDensity_timelocked{4}  = [];
SpikeTrials_timelocked{4}   = [];
SpikeTrials_windowed{4}     = [];
SpikeStats_windowed{4}      = [];
LFP{4}                      = [];
intervals{4}                = [];

% analyses
for ipatient = [1, 2, 3, 4]
    
%     % load settings
    config = pnh_setparams;
%     
    % read muse markers
    MuseStruct{ipatient} = readMuseMarkers(config{ipatient}, false);
    
    % align markers
    MuseStruct{ipatient} = alignMuseMarkersXcorr(config{ipatient}, MuseStruct{ipatient}, false);

    % add seizures 
    MuseStruct{ipatient} = updateMarkers(config{ipatient}, MuseStruct{ipatient}, {'CriseStart', 'CriseEnd'});
%            
%     % intervals between IEDs
%     intervals{ipatient} = inter_trial_intervals(config{ipatient}, MuseStruct{ipatient}, true);
% 
%     % read hypnogram for combining with windowed data (neccecary for stats but not used in article)
%     hypnogramMuseStats(config{ipatient}, MuseStruct{ipatient}, false);
% 
%     % read LFP data (trial-by-trial)
%     LFP{ipatient} = readLFP(config{ipatient}, MuseStruct{ipatient}, false);
% 
%     % Time-frequency of IEDs
%     TFR{ipatient} = TFRtrials(config{ipatient}, false);
%    
%     % add any artefacts added later
%     [MuseStruct{ipatient}] = updateBadMuseMarkers(config{ipatient}, MuseStruct{ipatient});
%     
%     % write parameters for spyking circus
%     writeSpykingCircusDeadfiles(config{ipatient}, MuseStruct_aligned{ipatient}, true);
%     writeSpykingCircusParameters(config{ipatient});
%     [filelist, sampleinfo, timestamps, hdr] = writeSpykingCircusFileList(config{ipatient}, false);

    % read spike data from Phy as one continuous trial
    SpikeRaw{ipatient} = readSpikeRaw_Phy(config{ipatient}, false);
     
    % segment into trials based on IED markers
    SpikeTrials_timelocked{ipatient}    = readSpikeTrials_MuseMarkers(config{ipatient}, MuseStruct{ipatient}, SpikeRaw{ipatient}, false);
 
    % spike density and stats vs. baseline
    SpikeDensity_timelocked{ipatient}   = spikeTrialDensity(config{ipatient}, SpikeTrials_timelocked{ipatient}, true);

    % create sliding timewindows
    config{ipatient}.spikewin.windowsize = 10; %Patient 1; 10: (0/total) = 7835/14826; 60: 352/2460
    for ipart = 1 : size(config{ipatient}.directorylist, 2)
        for idir = 1 : size(config{ipatient}.directorylist{ipart}, 2)
            temp = dir(fullfile(config{ipatient}.rawdir, config{ipatient}.directorylist{ipart}{idir}, ['*', config{ipatient}.LFP.channel{1}, '.ncs']));
            hdr  = ft_read_header(fullfile(config{ipatient}.rawdir, config{ipatient}.directorylist{ipart}{idir}, temp.name));
            MuseStruct{ipatient}{ipart}{idir}.markers.window__START__.synctime = 0 : (config{ipatient}.spikewin.windowsize - config{ipatient}.spikewin.windowsize * config{ipatient}.spikewin.windowoverlap) : hdr.nSamples/hdr.Fs - (config{ipatient}.spikewin.windowsize);
            MuseStruct{ipatient}{ipart}{idir}.markers.window__START__.clock    = seconds(MuseStruct{ipatient}{ipart}{idir}.markers.window__START__.synctime) + MuseStruct{ipatient}{ipart}{idir}.starttime;
            MuseStruct{ipatient}{ipart}{idir}.markers.window__START__.events   = size(MuseStruct{ipatient}{ipart}{idir}.markers.window__START__.synctime, 2);
            MuseStruct{ipatient}{ipart}{idir}.markers.window__END__.synctime   = MuseStruct{ipatient}{ipart}{idir}.markers.window__START__.synctime + config{ipatient}.spikewin.windowsize;
            MuseStruct{ipatient}{ipart}{idir}.markers.window__END__.clock      = MuseStruct{ipatient}{ipart}{idir}.markers.window__START__.clock + seconds(config{ipatient}.spikewin.windowsize);
            MuseStruct{ipatient}{ipart}{idir}.markers.window__END__.events     = size(MuseStruct{ipatient}{ipart}{idir}.markers.window__END__.synctime, 2);  
        end
    end

    config{ipatient}.muse.startmarker.window    = 'window__START__';
    config{ipatient}.muse.endmarker.window      = 'window__END__';
    config{ipatient}.spike.toi.window           = [0 0];
    config{ipatient}.spike.pad.window           = 0;
    config{ipatient}.spike.name                 = "window";
    config{ipatient}.spike.postfix              = '-windowed';
    SpikeTrials_windowed{ipatient}              = readSpikeTrials_MuseMarkers(config{ipatient}, MuseStruct{ipatient}, SpikeRaw{ipatient}, false);
    SpikeStats_windowed{ipatient}               = spikeTrialStats(config{ipatient}, SpikeTrials_windowed{ipatient}, true);

    % plot LFP timecourse examples for article
    config{ipatient}.plot.name      = {'patternA', 'patternB', 'patternC'};
    config{ipatient}.plot.ncols     = 4;
    config{ipatient}.plot.postfix   = "_IEDs";
    plotTimeCoursesExamples(config{ipatient});
    
    % plot seizure timecourse examples for article
    config{ipatient}.plot.name      = "seizureA";
    config{ipatient}.plot.ncols     = 1;    
    config{ipatient}.plot.postfix   = "_seizureA";
    plotTimeCoursesExamples(config{ipatient});
    
    % rasterplot & TFR (not in article, see Figure2.m & Figure3.m)
    config{ipatient}.plot.unit{1}   = ones(size(config{ipatient}.plot.unit{1})) * 0; % all in the same plot
    config{ipatient}.plot.name      = "SEIZURE";
    config{ipatient}.plot.ncols     = 1;    
    config{ipatient}.plot.postfix   = "_SEIZURE";
    plot_patterns_multilevel_examples(config{ipatient});
    
    % rasterplot & TFR (not in article, see Figure2.m & Figure3.m)
    config{ipatient}.plot.unit{1}   = ones(size(config{ipatient}.plot.unit{1})) * -1; % all individually
    if ipatient == 3
        config{ipatient}.plot.name   = {'FA','ES'};
    else
        config{ipatient}.plot.name   = {'PSW','FA','ES'};
    end
    config{ipatient}.plot.ncols     = 4;
    config{ipatient}.plot.postfix   = "_IEDs";
    plot_patterns_multilevel_examples2(config{ipatient});

    % plot firing rate change for each pattern (not in article, see Figure2.m & Figure3.m)
    plotstats(config{ipatient});
    
    % extract waveforms from raw data
    SpikeWaveforms = readSpikeWaveforms(config{ipatient}, SpikeRaw{ipatient}, false);

end

% summary of numbers and mode ISI of each pattern
t = table;

for ipatient = 1 : 4
    for markername = ["PSW", "FA", "ES"]
        if isfield(intervals{ipatient}.table, markername)
            fprintf('%s, Nodule %d: n=%d, λ=%0.0fms\n', markername, ipatient, size(intervals{ipatient}.table.(markername), 1), intervals{ipatient}.mode.(markername)*1000);
            
            fn = sprintf('%s_n', markername);
            t.(fn)(ipatient) = size(intervals{ipatient}.table.(markername), 1);
            
            fn = sprintf('%s_λ (ms)', markername);
            t.(fn)(ipatient) = intervals{ipatient}.mode.(markername)*1000;
           
        end
    end
end

% analyse periodicity c.f. Hirsch
analyze_periodicity(config, MuseStruct);
    
% Figure 2 for article
cfg_fig = [];
cfg_fig.colorscale = 0.2;
Figure2(cfg_fig);

% Figure 3 for article
Figure3;

% Table 2 for article
T = Table2(config, true);

% plot correlations between micro and macro

close all
fig = figure; hold;
mloc = [3, 2, 1; 3, 2, 1; 2, 1, 0];
xi = 0;
c = [];
labels = {};
for ipatient =  1 : 3
    for imarker = 1 : size(stat_TFR_micro{ipatient},2)
        xi = xi + 1;
        
        for icontact = 1 : length(stat_TFR_micro{ipatient}{imarker}.corrs.avg)
            xpos = (xi-1) * 3;
            ypos = (icontact-1) * 3;
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

clear stat_bl SpikeRateStats SpikeRateStats_bar SpikeRaw SpikeTrials tbl unit w_template w_data pt_template pt_data
tbl = table;
iunit = 1;

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
        tbl.nodule(iunit)   = ipatient;
        total               = length(SpikeRateStats{ipatient}.isi{itemp});
        short               = sum(SpikeRateStats{ipatient}.isi{itemp} <= 2);
        long                = sum(SpikeRateStats{ipatient}.isi{itemp} > 2);
        tbl.unit(iunit)     = itemp;
        tbl.total(iunit)    = total;
        
        if short/total * 100 < 1
            tbl.percRPV{iunit} = sprintf('$%.2f^s$',round(short/total * 100,2));
        else
            tbl.percRPV{iunit} = sprintf('$%.2f^m$',round(short/total * 100,2));
        end
        
        for ilabel = 1 : size(SpikeRateStats_bar{ipatient}.clusterstat,2)
            
            if isfield(SpikeRateStats_bar{ipatient}.clusterstat{ilabel}{itemp},'posclusters')
                if ~isempty(SpikeRateStats_bar{ipatient}.clusterstat{ilabel}{itemp}.posclusters)
                    if SpikeRateStats_bar{ipatient}.clusterstat{ilabel}{itemp}.posclusters(1).prob < 0.025 % stats are sorted to most significant
                        if SpikeRateStats_bar{ipatient}.clusterstat{ilabel}{itemp}.posclusters(1).prob < 0.001
                            tbl.([config{ipatient}.name{ilabel},'_increase']){iunit} = sprintf('$%.0f^{xxx}$',SpikeRateStats_bar{ipatient}.clusterstat{ilabel}{itemp}.maxcluster.perc{1});
                        elseif SpikeRateStats_bar{ipatient}.clusterstat{ilabel}{itemp}.posclusters(1).prob < 0.01
                            tbl.([config{ipatient}.name{ilabel},'_increase']){iunit} = sprintf('$%.0f^{xx}$', SpikeRateStats_bar{ipatient}.clusterstat{ilabel}{itemp}.maxcluster.perc{1});
                        else
                            tbl.([config{ipatient}.name{ilabel},'_increase']){iunit} = sprintf('$%.0f^{x}$',  SpikeRateStats_bar{ipatient}.clusterstat{ilabel}{itemp}.maxcluster.perc{1});
                        end
                    else
                        tbl.([config{ipatient}.name{ilabel},'_increase']){iunit} = 'n.s.';
                    end
                    
                else
                    tbl.([config{ipatient}.name{ilabel},'_increase']){iunit} = 'n.s.';
                end
            else
                tbl.([config{ipatient}.name{ilabel},'_increase']){iunit} = 'n.s.';
            end
            
            if isfield(SpikeRateStats_bar{ipatient}.clusterstat{ilabel}{itemp},'negclusters')
                if ~isempty(SpikeRateStats_bar{ipatient}.clusterstat{ilabel}{itemp}.negclusters)
                    if SpikeRateStats_bar{ipatient}.clusterstat{ilabel}{itemp}.negclusters(1).prob < 0.025 % stats are sorted to most significant
                        if SpikeRateStats_bar{ipatient}.clusterstat{ilabel}{itemp}.negclusters(1).prob < 0.001
                            tbl.([config{ipatient}.name{ilabel},'_decrease']){iunit} = sprintf('$%.0f^{xxx}$',SpikeRateStats_bar{ipatient}.clusterstat{ilabel}{itemp}.mincluster.perc{1});
                        elseif SpikeRateStats_bar{ipatient}.clusterstat{ilabel}{itemp}.negclusters(1).prob < 0.01
                            tbl.([config{ipatient}.name{ilabel},'_decrease']){iunit} = sprintf('$%.0f^{xx}$',SpikeRateStats_bar{ipatient}.clusterstat{ilabel}{itemp}.mincluster.perc{1});
                        else
                            tbl.([config{ipatient}.name{ilabel},'_decrease']){iunit} = sprintf('$%.0f^{x}$',SpikeRateStats_bar{ipatient}.clusterstat{ilabel}{itemp}.mincluster.perc{1});
                        end
                    else
                        tbl.([config{ipatient}.name{ilabel},'_decrease']){iunit} = 'n.s.';
                    end
                else
                    tbl.([config{ipatient}.name{ilabel},'_decrease']){iunit} = 'n.s.';
                end
            else
                tbl.([config{ipatient}.name{ilabel},'_decrease']){iunit} = 'n.s.';
            end
            
            
            if corrs{ipatient}{ilabel}{itemp}.p < 0.001
                tbl.([config{ipatient}.name{ilabel},'_corr']){iunit}  = sprintf('$%.2f^{xxx}$',corrs{ipatient}{ilabel}{itemp}.rho);
            elseif corrs{ipatient}{ilabel}{itemp}.p < 0.01
                tbl.([config{ipatient}.name{ilabel},'_corr']){iunit}  = sprintf('$%.2f^{xx}$',corrs{ipatient}{ilabel}{itemp}.rho);
            elseif corrs{ipatient}{ilabel}{itemp}.p < 0.05
                tbl.([config{ipatient}.name{ilabel},'_corr']){iunit}  = sprintf('$%.2f^{x}$',corrs{ipatient}{ilabel}{itemp}.rho);
            else
                tbl.([config{ipatient}.name{ilabel},'_corr']){iunit}  = 'n.s';
            end
            
            %             tbl.([config{ipatient}.name{ilabel},'_corr']) =
            
        end
        
        % zero crossing function
        zci                 = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);
        
        % take the first file to extract the header of the data
        hdr_fname           = fullfile(temp(1).folder,temp(1).name);
        hdr                 = ft_read_header(hdr_fname);
        
        % read muse markers
        [MuseStruct_micro, MuseStruct_macro]    = readMuseMarkers(config{ipatient}, false);
        
        % align Muse markers according to peaks and detect whether they contain artefacts
        [MuseStruct_micro, MuseStruct_macro]    = alignMuseMarkers(config{ipatient},MuseStruct_micro, MuseStruct_macro, false);
        
        % read raw spike data from SC, and segment into trials
        [SpikeRaw, SpikeTrials]                 = readSpykingCircus(config{ipatient}, MuseStruct_micro, false);
        
        [spike_LFP]  = spikeLFP(config{ipatient}, SpikeRaw, false);
        
        % define width of spike to discard overlapping spike times
        w_pre  = config{ipatient}.spike.pre  * spike_LFP{itemp}.fsample;
        w_post = config{ipatient}.spike.post * spike_LFP{itemp}.fsample;
        
        % select non-overlapping spikes and average
        spike_LFP_sel = spike_LFP{itemp};
        i = 1;
        for trialnr = 1 : size(spike_LFP{itemp}.trial,2)
            if size(find(SpikeRaw.sample{itemp} < SpikeRaw.sample{itemp}(trialnr)+w_post & SpikeRaw.sample{itemp} > SpikeRaw.sample{itemp}(trialnr)-w_pre),1) == 1
                fprintf('Adding unit %d of %d: trial %d of %d\n',itemp,size(spike_LFP,2),trialnr,size(spike_LFP{itemp}.trial,2));
                spike_LFP_sel.trial{i} = spike_LFP{itemp}.trial{trialnr};
                i = i + 1;
            else
                fprintf('Skipping unit %d of %d: trial %d of %d because it overlaps with other spikes\n',itemp,size(spike_LFP,2),trialnr,size(spike_LFP{itemp}.trial,2));
            end
        end
        
        % average waveshapes and apply baseline correction
        spike_LFP_avg{itemp}    = ft_timelockanalysis([],spike_LFP_sel);
        %         cfgtemp                 = [];
        %         cfgtemp.demean          = 'yes';
        %         cfgtemp.baselinewindow  = [-0.001, -0.0005];
        %         spike_LFP{itemp}        = ft_preprocessing(cfgtemp,spike_LFP{itemp});
        
        % extract timecourse and timeaxis of waveshape
        tempsel                 = spike_LFP_avg{itemp}.avg;
        temptime                = spike_LFP_avg{itemp}.time;
        
        % interpolate template
        temptime_int            = linspace(temptime(1),temptime(end),10000);
        tempsel_int             = pchip(temptime,tempsel,temptime_int);
        
        % save data for figures later (ugly code);
        tempsel_data_forfig{iunit}   = tempsel;
        temptime_data_forfig{iunit}  = temptime;
        
        %% plot spike according to raw data
        
        % initialize figure
        fig = figure;
        
        % find positive peak - i.e. spike peak
        [Ypos,Xpos]             = findpeaks(tempsel_int,temptime_int,'NPeaks',1,'SortStr','descend');
        
        % 0 based on data timeaxis, to select trough
        izero                   = find(temptime_int > 0,1,'first');
        
        % find negative peak after t=0, i.e. trough
        [Yneg,Xneg]             = findpeaks(-tempsel_int(izero:end),temptime_int(izero:end),'NPeaks',1,'SortStr','descend');
        
        % plot
        subplot(1,2,1); hold;
        plot(temptime_int,tempsel_int);
        plot([Xpos,Xneg(1)], [Yneg(1)*0.3, Yneg(1)*0.3],'-o','Color',[1 0 0],'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]);
        
        tbl.pt_data{iunit}      = sprintf('%.2f',abs(Xpos-Xneg)*1000);
        pt_data(iunit)          = abs(Xpos-Xneg)*1000;
        %         midline                     = min(tempsel_int) + (max(tempsel_int) - min(tempsel_int)) / 2 ;
        midline                     = max(tempsel_int) / 2 ;
        indx                        = zci(tempsel_int - midline);
        tbl.width_data{iunit}   = sprintf('%.2f',diff(temptime_int(indx))*1000);
        w_data(iunit)           = diff(temptime_int(indx))*1000;
        
        plot([temptime_int(indx(1)),temptime_int(indx(2))], [midline, midline],'-o','Color',[0 0 1],'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]);
        title('Data');
        
        
        %% plot template
        
        % load data
        temp                = dir(fullfile(config{ipatient}.datasavedir,[config{ipatient}.prefix,'all_data_',config{ipatient}.circus.channel{1}(1:end-2),'_*.ncs']));
        tempsel             = squeeze(SpikeRaw.template(itemp,SpikeRaw.template_maxchan(itemp),:));
        temptime            = ((1:size(SpikeRaw.template,3))/hdr.Fs*1000)';
        
        % interpolate template
        temptime_int        = linspace(temptime(1),temptime(end),10000);
        tempsel_int         = pchip(temptime,tempsel,temptime_int);
        
        % save data for figures later (ugly code);
        tempsel_template_forfig{iunit}   = tempsel;
        temptime_template_forfig{iunit}  = temptime;
        
        % find positive peak - i.e. spike peak
        [Ypos,Xpos] = findpeaks(tempsel_int,temptime_int,'NPeaks',1,'SortStr','descend');
        
        % t=0 based on template
        izero = find((temptime_int - temptime_int(end)/2) > 0,1,'first');
        
        % find negative peak after t=0, i.e. trough
        [Yneg,Xneg]         = findpeaks(-tempsel_int(izero:end),temptime_int(izero:end),'NPeaks',1,'SortStr','descend');
        
        subplot(1,2,2); hold;
        plot(temptime_int,tempsel_int);
        plot([Xpos,Xneg(1)], [Yneg(1)*0.3, Yneg(1)*0.3],'-o','Color',[1 0 0],'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]);
        
        tbl.pt_template{iunit}      = sprintf('%.2f',abs(Xpos-Xneg)*1000);
        pt_template(iunit)          = abs(Xpos-Xneg)*1000;
        %         midline                 = min(tempsel_int) + (max(tempsel_int) - min(tempsel_int)) / 2 ;
        midline                     = max(tempsel_int) / 2 ;
        indx                        = zci(tempsel_int - midline);
        tbl.width_template{iunit}   = sprintf('%.2f',diff(temptime_int(indx))*1000);
        w_template(iunit)           = diff(temptime_int(indx))*1000;
        
        plot([temptime_int(indx(1)),temptime_int(indx(2))], [midline, midline],'-o','Color',[0 0 1],'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]);
        title('Template');
        axis tight
        
        % save figure
        fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
        set(fig,'PaperOrientation','landscape');
        set(fig,'PaperUnits','normalized');
        set(fig,'PaperPosition', [0 0 1 1]);
        print(fig, '-dpdf', fullfile(config{ipatient}.imagesavedir,[config{ipatient}.prefix,'spike_morphology_U',num2str(itemp),'.pdf']),'-r600');
        
        % keep counting units across patient for table
        iunit = iunit + 1;
    end
end

tbl.displayed{tbl.nodule == 1 & tbl.unit == 2} = 'A';
tbl.displayed{tbl.nodule == 2 & tbl.unit == 8} = 'B';
tbl.displayed{tbl.nodule == 2 & tbl.unit == 5} = 'C';
tbl.displayed{tbl.nodule == 3 & tbl.unit == 7} = 'D';

% tbl = sortrows(tbl,{'nodule','SUA'},{'ascend','descend'});

% writetable(tbl,'/network/lustre/iss01/charpier/stephen.whitmarsh/WhitmarshEpilepsy/unittables.xls')
writetable(tbl,'/network/lustre/iss01/charpier/stephen.whitmarsh/WhitmarshEpilepsy/unittables_spikewaveshapeaverage.xls')


%% Make figure like Bartho 2000

fig = figure;

subplot(1,2,1); hold;
i = strcmp(tbl.displayed, 'A') | strcmp(tbl.displayed, 'B') | strcmp(tbl.displayed, 'C') | strcmp(tbl.displayed, 'D');
scatter(w_data(i)/1000,pt_data(i)/1000,120,'r','d','filled');
xlabel('Half-amplitude duration');
ylabel('Trough to peak time');

% i = strcmp(tbl.SUA, 'MUA');
scatter(w_data/1000,pt_data/1000,60,'k','linewidth',2);
%
% i = strcmp(tbl.SUA, 'SUA');
% scatter(w(i)/1000,pt(i)/1000,60,'k','filled');

for i = 1 : size(w_data,2)
    text(w_data(i)/1000+0.00001,pt_data(i)/1000+0.00005,sprintf('%d-%d',tbl.nodule(i),tbl.unit(i)));
end

title('Data');
axis square

subplot(1,2,2); hold;
i = strcmp(tbl.displayed, 'A') | strcmp(tbl.displayed, 'B') | strcmp(tbl.displayed, 'C') | strcmp(tbl.displayed, 'D');
scatter(w_template(i)/1000,pt_template(i)/1000,120,'r','d','filled');
xlabel('Half-amplitude duration');
ylabel('Trough to peak time');

% i = strcmp(tbl.SUA, 'MUA');
scatter(w_template/1000,pt_template/1000,60,'k','linewidth',2);
%
% i = strcmp(tbl.SUA, 'SUA');
% scatter(w(i)/1000,pt(i)/1000,60,'k','filled');

for i = 1 : size(w_template,2)
    text(w_template(i)/1000+0.003,pt_template(i)/1000+0.015,sprintf('%d-%d',tbl.nodule(i),tbl.unit(i)));
end
title('Template');
axis square

% print ISI to file
fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
set(fig,'PaperOrientation','landscape');
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition', [0 0 1 1]);
print(fig, '-dpdf', '/network/lustre/iss01/charpier/stephen.whitmarsh/WhitmarshEpilepsy/scatter_spikes_time_template.pdf','-r600');

%% separate units according to trough-peak

inti = find(pt_template <= 500);
pyri = find(pt_template > 500);

fig = figure;
subplot(2,2,1); hold;
for i = inti
    x = temptime_data_forfig{i};
    y = tempsel_data_forfig{i};
    y = y / max(y);
    plot(x, y,'k');
end
title('Data Interneurons');
axis tight;

subplot(2,2,2); hold;
for i = inti
    x = temptime_template_forfig{i};
    y = tempsel_template_forfig{i};
    y = y / max(y);
    plot(x, y,'k');
end
title('Template Interneurons');
axis tight;

subplot(2,2,3); hold;
for i = pyri
    x = temptime_data_forfig{i};
    y = tempsel_data_forfig{i};
    y = y / max(y);
    plot(x, y,'k');
end
title('Data Pyramidal');
axis tight;

subplot(2,2,4); hold;
for i = pyri
    x = temptime_template_forfig{i};
    y = tempsel_template_forfig{i};
    y = y / max(y);
    plot(x, y,'k');
end
title('Template Pyramidal');
axis tight;

% print to file
fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
set(fig,'PaperOrientation','landscape');
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition', [0 0 1 1]);
print(fig, '-dpdf', '/network/lustre/iss01/charpier/stephen.whitmarsh/WhitmarshEpilepsy/fig4_average_waveforms.pdf','-r600');


%% plot example of individual waveforms of single unit

% 2-8 (interneuron) and 1-2 (pyramidal)

ipatient = 2;
itemp = 8;

ipatient = 1;
itemp = 2;

for ipatient = [1,2]
    if ipatient == 1
        itemp = 2;
    end
    if ipatient == 2
        itemp = 8;
    end
    
    config = pnh_setparams([]);
    
    % read muse markers
    [MuseStruct_micro, MuseStruct_macro]    = readMuseMarkers(config{ipatient}, false);
    
    % align Muse markers according to peaks and detect whether they contain artefacts
    [MuseStruct_micro, MuseStruct_macro]    = alignMuseMarkers(config{ipatient},MuseStruct_micro, MuseStruct_macro, false);
    
    % read raw spike data from SC, and segment into trials
    [SpikeRaw, ~]                           = readSpykingCircus(config{ipatient}, MuseStruct_micro, false);
    
    [spike_LFP]                             = spikeLFP(config{ipatient}, SpikeRaw, false);
    
    cfgtemp                 = [];
    cfgtemp.demean          = 'yes';
    cfgtemp.baselinewindow  = 'all'; % baseline to peak to normalize peak-to-peak
    spike_LFP{itemp}        = ft_preprocessing(cfgtemp,spike_LFP{itemp});
    
    % define width of spike to discard overlapping spike times
    w_pre  = config{ipatient}.spike.pre  * spike_LFP{itemp}.fsample;
    w_post = config{ipatient}.spike.post * spike_LFP{itemp}.fsample;
    
    % select non-overlapping spikes and average
    clear spike
    i = 1;
    for trialnr = 1 : size(spike_LFP{itemp}.trial,2)
        if size(find(SpikeRaw.sample{itemp} < SpikeRaw.sample{itemp}(trialnr)+w_post & SpikeRaw.sample{itemp} > SpikeRaw.sample{itemp}(trialnr)-w_pre),1) == 1
            fprintf('Adding unit %d of %d: trial %d of %d\n',itemp,size(spike_LFP,2),trialnr,size(spike_LFP{itemp}.trial,2));
            spike(i,:) = spike_LFP{itemp}.trial{trialnr};
            i = i + 1;
        else
            fprintf('Skipping unit %d of %d: trial %d of %d because it overlaps with other spikes\n',itemp,size(spike_LFP,2),trialnr,size(spike_LFP{itemp}.trial,2));
        end
    end
    
    i = var(spike,0,2) < std(var(spike,0,2))*0.5;
    spike_sel = spike(i,:);
    
    r = randi(size(spike_sel,1),[1000,1]);
    % r = randi(size(spike_sel,1),[size(spike_sel,1),1]);
    
    fig = figure; hold;
    for i = r
        plot(spike_LFP{itemp}.time{1},spike_sel(i,:),'k');
    end
    plot(spike_LFP{itemp}.time{1},mean(spike),'r','linewidth',3);
    
    
    % print to file
    fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
    set(fig,'PaperOrientation','landscape');
    set(fig,'PaperUnits','normalized');
    set(fig,'PaperPosition', [0 0 1 1]);
    print(fig, '-dpdf', ['/network/lustre/iss01/charpier/stephen.whitmarsh/WhitmarshEpilepsy/example_waveforms_nodule',num2str(ipatient),'_unit',num2str(itemp),'.pdf'],'-r600');
end

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
        isi = diff(SpikeRaw{ipatient}.sample{itemp}) / hdr.Fs * 1000;
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

