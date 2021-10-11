%% Analysis script for Valerio Frazzini's study on electrophysiology of pathoanatomy
%
% (c) Stephen Whitmarsh, stephen.whitmarsh@gmail.com
%
% requires bandpassFilter.m from Mario
% requires releaseDec2015 from Neuralynx website

%% Add path

restoredefaultpath
if isunix
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/git/fieldtrip
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/projects/pnh/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/shared/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/shared/utilities
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/subaxis    
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/fieldtrip    
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/sigstar-master
    addpath(genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/scripts/releaseDec2015/'));    
    addpath(genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/scripts/epishare-master'));
    addpath(genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/scripts/SPIKY_apr_2021'))    
end

if ispc
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\git\fieldtrip
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\projects\pnh
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\shared
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\shared\utilities
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external\altmany-export_fig-8b0ba13\
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external\subaxis    
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external\fieldtrip    
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external\sigstar-master
    addpath(genpath('\\lexport\iss01.charpier\analyses\stephen.whitmarsh\scripts\epishare-master'));
    addpath(genpath('\\lexport\iss01.charpier\analyses\stephen.whitmarsh\scripts\SPIKY_apr_2021'));
    addpath          \\lexport\iss01.charpier\analyses\stephen.whitmarsh\scripts\MatlabImportExport_v6.0.0
end
ft_defaults

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
for ipatient = 1 : 4
    
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
    % read LFP data (trial-by-trial)
    config{ipatient}.LFP.name = "SEIZURE";
    LFP{ipatient} = readLFP(config{ipatient}, MuseStruct{ipatient}, true);
    
    
% 
%     % Time-frequency of IEDs
%     TFR{ipatient} = TFRtrials(config{ipatient}, false);
%    
%     % add any artefacts added later
%     [MuseStruct{ipatient}] = updateBadMuseMarkers(config{ipatient}, MuseStruct{ipatient});
%     
%     % write parameters for spyking circus
    writeSpykingCircusDeadfiles(config{ipatient}, MuseStruct_aligned{ipatient}, true);
    writeSpykingCircusParameters(config{ipatient});
    [filelist, sampleinfo, timestamps, hdr] = writeSpykingCircusFileList(config{ipatient}, false);

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

    % plot LFP timecourse examples for Figure 1
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
    plot_patterns_multilevel_examples(config{ipatient});

    % plot firing rate change for each pattern (not in article, see Figure2.m & Figure3.m)
    plotstats(config{ipatient});
    
    % extract waveforms from raw data
    if ipatient == 1
        SpikeRaw{ipatient}{1}.template_maxchan(2) = 2;
        SpikeWaveforms = readSpikeWaveforms(config{ipatient}, SpikeRaw{ipatient}, false);
        for itrial = 1 : size(SpikeWaveforms{1}{2}.trial, 2)
            SpikeWaveforms{1}{2}.trial{itrial} = -SpikeWaveforms{1}{2}.trial{itrial};
        end
        fname = fullfile(config{ipatient}.datasavedir, [config{ipatient}.prefix, 'spike_waveform.mat']);
        save(fname, 'SpikeWaveforms', '-v7.3');
    else
        SpikeWaveforms = readSpikeWaveforms(config{ipatient}, SpikeRaw{ipatient}, false);
    end

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
    
% plot overview pictures
config                  = pnh_setparams;
config{1}.plot.name     = 'ES';
config{1}.plot.macro    = {'_1pNs_1', '_1pNs_2', '_1pNs_3'};
config{1}.plot.micro    = {'m1pNs_4', 'm1pNs_6'};
config{1}.plot.ipart    = 1;
config{1}.plot.postfix  = [];
config{1}.plot.ievent   = 14;

config{2}.plot.name     = 'ES';
config{2}.plot.macro    = {'_Casd_1', '_Casd_2', '_Casd_3'};
config{2}.plot.micro    = {'mCasd_1', 'mCasd_2'};
config{2}.plot.ipart    = 1;
config{2}.plot.postfix  = [];
config{2}.plot.ievent   = 14;

config{3}.plot.name     = 'ES';
config{3}.plot.macro    = {'_TNmi_1', '_TNmi_2', '_TNmi_3'};
config{3}.plot.micro    = {'mTNmi_2', 'mTNmi_3'};
config{3}.plot.ipart    = 1;
config{3}.plot.postfix  = [];
config{3}.plot.ievent   = 1;

config{4}.plot.name     = 'ES';
config{4}.plot.macro    = {'_LMI1_1', '_LMI1_2', '_LMI1_3'};
config{4}.plot.micro    = {'mLMI1_1', 'mLMI1_7'};
config{4}.plot.ipart    = 1;
config{4}.plot.postfix  = [];
config{4}.plot.ievent   = 14;

config                  = pnh_setparams;
config{5} = config{1};
config{5}.plot.name     = 'SEIZURE';
config{5}.plot.macro    = {'_1pNs_1', '_1pNs_2', '_1pNs_3'};
config{5}.plot.micro    = {'m1pNs_4', 'm1pNs_6'};
config{5}.plot.ipart    = 1;
config{5}.plot.postfix  = [];
config{5}.plot.ievent   = 1;
Figure0(config{5}, MuseStruct{1});

config                  = pnh_setparams;
config{6} = config{4};
config{6}.plot.name     = 'SEIZURE';
config{6}.plot.macro    = {'_LMI1_1', '_LMI1_2', '_LMI1_3'};
config{6}.plot.micro    = {'mLMI1_1', 'mLMI1_7'};
config{6}.plot.ipart    = 1;
config{6}.plot.postfix  = [];
config{6}.plot.ievent   = 2;
Figure0(config{6}, MuseStruct{4});


for ipatient = 1 : 4
    for ievent = 10 : 10 : 100
        
        config{ipatient}.plot.ievent = ievent;
        % read muse markers
        MuseStruct{ipatient} = readMuseMarkers(config{ipatient}, false);
        
        % align markers
        MuseStruct{ipatient} = alignMuseMarkersXcorr(config{ipatient}, MuseStruct{ipatient}, false);
        
        Figure0(config{ipatient}, MuseStruct{ipatient});
        close all
        
    end
end

% Figure 2 for article
cfg_fig = [];
cfg_fig.colorscale = 0.2;
Figure2(cfg_fig);

% Figure 3 for article
Figure3;

% Figure 4 for article
Figure4;

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

