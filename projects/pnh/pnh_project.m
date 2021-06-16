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
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/projects/hspike/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/shared/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/shared/utilities
    addpath(genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/scripts/releaseDec2015/'));
    addpath(genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/epishare-master'));
end

if ispc
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\fieldtrip
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\projects\hspike
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\shared
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\shared\utilities
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external\altmany-export_fig-8b0ba13\
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\MatlabImportExport_v6.0.0
    addpath(genpath('\\lexport\iss01.charpier\analyses\stephen.whitmarsh\epishare-master'));
end
ft_defaults

feature('DefaultCharacterSet', 'CP1252') % To fix bug for weird character problems in reading neurlynx

% load settings
config = pnh_setparams;


%% General analyses


for ipatient =  1 : 3
    
    
    % align Muse markers according to peaks and detect whether they contain artefacts
    [MuseStruct_micro, MuseStruct_macro]    = alignMuseMarkers(config{ipatient},MuseStruct_micro, MuseStruct_macro, false);
%     [MuseStruct_micro, MuseStruct_macro]    = MuseMarkers_update_filepath(config{ipatient},MuseStruct_micro, MuseStruct_macro);
    
%     % read LFP data
    [dat_micro, dat_macro] = readLFP(config{ipatient},MuseStruct_micro, MuseStruct_macro, false);
%     
%     % plot LFP timecourse examples for article
%     plotTimeCourses(config{ipatient});
%     )
%     % plot LFP data
%     [FFT_micro_trials,TFR_micro_trials,TFR_macro_trials,stat_TFR_micro] = plotLFP(config{ipatient}, dat_micro, dat_macro, true);
%     
    % write data concatinated for SC, and update config with sampleinfo
    config{ipatient}                        = writeSpykingCircus(config{ipatient}, MuseStruct_micro, false);
    
    % read raw spike data from SC, and segment into trials
    [SpikeRaw, SpikeTrials]                 = readSpykingCircus(config{ipatient}, MuseStruct_micro, false);
    
    % read and plot spikerate overview, and get the stats
    [SpikeRateStats{ipatient}, stats_bar{ipatient}, sdf_orig_out{ipatient}, sdf_bar_out] = spikeratestats(config{ipatient}, SpikeRaw, SpikeTrials, false);
    
    % read and plot LFP of spike events
%     [spike_LFP]  = spikeLFP(config{ipatient},SpikeRaw, false);
    
end


%% inter-event timings
fig = figure;
subplots  = [3, 2, 1, 6, 5, 4, 8, 7];
binwidths = [1, 1, 0.25, 1, 1, 0.25, 1, 0.25];
binlimits = [60, 60, 15, 60, 60, 15, 60, 15];
i = 1;

for ipatient =  1 : 3
    
    config = pnh_setparams([]);
    
    % read muse markers
    [MuseStruct_micro, MuseStruct_macro]    = readMuseMarkers(config{ipatient}, false);
    
    % align Muse markers according to peaks and detect whether they contain artefacts
    [MuseStruct_micro, MuseStruct_macro]    = alignMuseMarkers(config{ipatient},MuseStruct_micro, MuseStruct_macro, false);
    
    
    % read LFP data
    [intervals] = inter_trial_intervals(config{ipatient},MuseStruct_micro,true);
    s = fieldnames(intervals);
    
    for ss = s'
        subplot(3,3,subplots(i));
        histogram(intervals.(cell2mat(ss))(:,4),'BinLimits',[0 binlimits(i)],'BinWidth',binwidths(i),'EdgeColor','black','facecolor','k');
        me = nanmedian(intervals.(cell2mat(ss))(:,4));
        mo = mode(intervals.(cell2mat(ss))(:,4));
        m  = nanmean(intervals.(cell2mat(ss))(:,4));
        sd = nanstd(intervals.(cell2mat(ss))(:,4));
        axis tight
%         title(sprintf('Nodule: %d, Pattern: %s, Median: %1.2f, Mode: %1.2f, Mean: %1.2f, SD: %1.2f',ipatient,cell2mat(ss),me,mo,m,sd));
        title(sprintf('Nodule: %d, Pattern: %s',ipatient,cell2mat(ss)));
        xlabel('Seconds');
        ylabel('Count');
        box off
        
        %         % print ISI to file
        %         fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
        %         set(fig,'PaperOrientation','landscape');
        %         set(fig,'PaperUnits','normalized');
        %         set(fig,'PaperPosition', [0 0 1 1]);
        %         xlabel('Seconds');
        %         ylabel('Count');
        %         box off
        %
        %         print(fig, '-dpdf', fullfile(config{ipatient}.imagesavedir,['N',num2str(ipatient),'_inter-trial-intervals_',cell2mat(ss)]));
        i = i + 1;
    end
    
end
% print ISI to file
fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
set(fig,'PaperOrientation','landscape');
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition', [0 0 1 1]);

print(fig, '-dpdf', fullfile(config{ipatient}.imagesavedir,'Inter-trial-intervals'));


%% revised calculation of unit firing charateristics
        
for ipatient = 1 : 3
    
    [SpikeRaw{ipatient},~]  = readSpykingCircus(config{ipatient}, MuseStruct_micro, false);
    
    % get samplerate
    temp                    = dir(fullfile(config{ipatient}.datasavedir,[config{ipatient}.prefix,'all_data_',config{ipatient}.circus.channel{1},'.ncs']));
    hdr_fname               = fullfile(temp(1).folder,temp(1).name);
    hdr                     = ft_read_header(hdr_fname); % take the first file to extract the header of the data
    
    % redefine trials to 60-second windows for ISI
    window = 60;
    cfgtemp                 = [];
    cfgtemp.trl             = (1 : hdr.Fs*window : hdr.nSamples)';
    cfgtemp.trl(:,2)        = cfgtemp.trl(:,1) + hdr.Fs*window;
    cfgtemp.trl(:,3)        = zeros(size(cfgtemp.trl,1),1);
    %         cfgtemp.trl             = cfgtemp.trl(1:end-1,:);
    %         cfgtemp.trl             = cfgtemp.trl(randi(size(cfgtemp.trl,1),1000,1),:);    % select a subsection, to reduce memoryload
    cfgtemp.trlunit         = 'samples';
    cfgtemp.hdr             = hdr;
    SpikeTrials{ipatient}   = ft_spike_maketrials(cfgtemp,SpikeRaw{ipatient});
end

% ISI descriptives
clear stats
for ipatient = 1 : 3    
    
     for itemp = 1 : length(SpikeTrials{ipatient}.label)
        
        clear trialavg_isi
        isi_intraburst    = [];
        isi_interburst    = [];
        i           = 1;
        isi_pooled  = [];
        
        for itrial = unique(SpikeTrials{ipatient}.trial{itemp})
            % get timings and ISIs per trial
            indx            = SpikeTrials{ipatient}.trial{itemp} == itrial;
            t               = SpikeTrials{ipatient}.time{itemp}(indx);
            isi_all         = diff(t);
            
            % counting bursts as in Colder et al. 1996, & Staba et al. 2002
            indx            = isi_all < 0.01; % ask Stephane
            burstindx       = zeros(size(indx));
            toremove        = [];
            
            for burstlength = 1 : 10
                
                pattern     = [false, true(1,burstlength), false];
                bindx       = strfind(indx, pattern);
                
                if ~isempty(bindx)
                    burstindx(bindx+1) = burstlength; % note +1 because pattern starts with zero
                    fprintf('Found %d bursts of length %d in trial %d \n',length(bindx), burstlength, itrial);
                    
                    % add to list to correct for bursts
                    for ii = 1 : size(bindx,2)
                        
                        % remove all but first spike (at +1)
                        toremove = [toremove, bindx(ii)+2:bindx(ii)+2+burstlength-1]; % burstlength = 1; 0 1 0 -> 0 1 x 0
                        
                        % add ISI within bursts
                        isi_intraburst = [isi_intraburst, isi_all(bindx(ii)+1:bindx(ii)+burstlength-1)];
                        
                    end
                end
                
                stats{ipatient}.burstsum{itemp}(itrial, burstlength) = sum(length(bindx));
                
            end
            
            % concatinate ISIs, but only between bursts (not within)
            t_interburst               = t(burstindx ~= 0);
            isi_interburst             = [isi_interburst, diff(t_interburst)];
            
            % remove subsequenct APs after first AP of a burst
            t_corrected                 = t;
            t_corrected(toremove)       = [];
            isi_corrected               = diff(t_corrected);
            
            % basic descriptives
            trialavg_isi(i)             = nanmean(isi_all);
            trialfreq(i)                = 1/nanmean(isi_all);
            spikecount(i)               = size(t,2);
            spikecount_corrected(i)     = size(t_corrected,2);
            
            % according to Ponce-Alvarez, 2010
            x                           = isi_corrected(1:end-1) ./ isi_corrected(2:end);
            CV2_instant                 = 2 * abs(x - 1) ./ (x + 1);
            CV2_trial(i)                = mean(CV2_instant);
            
            x                           = isi_intraburst(1:end-1) ./ isi_intraburst(2:end);
            CV2_intraburst_instant      = 2 * abs(x - 1) ./ (x + 1);
            CV2_intraburst_trial(i)     = mean(CV2_intraburst_instant);
            
            LV_instant                  = 3 * (x - 1).^2 ./ (x + 1).^2;
            LV_trial(i)                 = mean(LV_instant);
            IR_instant                  = abs(log(x));
            IR_trial(i)                 = mean(IR_instant);
            SI_instant                  = 0.5 * log((x+1).^2/(4*x));
            SI_trial(i)                 = mean(SI_instant);
            
            % concatinate ISIS over trials, corrected for bursts: for pooled CV
            isi_pooled                  = [isi_pooled, isi_corrected];
            
            % calculate CV per trial for averged CV
            CV_trial(i)                 = nanstd(isi_corrected) / nanmean(isi_corrected);
            
            % short vs long ISIs for BI
            short(i)                    = sum(isi_all < 0.010);
            long(i)                     = sum(isi_all < 0.100);
            
            i = i + 1;
        end
        
        % get stats per sleep stage, over trials
        stats{ipatient}.isi_intraburst{itemp}            = isi_intraburst;
        stats{ipatient}.isi_interburst{itemp}            = isi_interburst;
        stats{ipatient}.burst_trialsum{itemp}            = sum(stats{ipatient}.burstsum{itemp});
        stats{ipatient}.mean_freq{itemp}                 = nanmean(trialfreq);
        stats{ipatient}.stdev_freq{itemp}                = nanstd(trialfreq);
        [N,EDGES]                                        = histcounts(trialfreq,'BinWidth',0.5);
        [M,I]                                            = max(N);
        stats{ipatient}.mode_freq{itemp}                 = mean(EDGES(I:I+1));
        stats{ipatient}.mean_isi{itemp}                  = nanmean(trialavg_isi);
        stats{ipatient}.var_isi{itemp}                   = nanstd(isi_pooled)^2;
        stats{ipatient}.burstindex{itemp}                = sum(short) / sum(long);
        stats{ipatient}.FF{itemp}                        = nanstd(spikecount_corrected)^2 / nanmean(spikecount_corrected);
        stats{ipatient}.spikecount{itemp}                = sum(spikecount);
        stats{ipatient}.spikecount_corrected{itemp}      = sum(spikecount_corrected);
        stats{ipatient}.CV_pooled{itemp}                 = nanstd(isi_pooled)   / nanmean(isi_pooled);
        stats{ipatient}.CV_trialavg{itemp}               = nanmean(CV_trial);
        stats{ipatient}.CV2_trialavg{itemp}              = nanmean(CV2_trial);
        stats{ipatient}.CV2_intraburst_trialavg{itemp}   = nanmean(CV2_intraburst_trial);
        stats{ipatient}.LV_trialavg{itemp}               = nanmean(LV_trial);
        stats{ipatient}.IR_trialavg{itemp}               = nanmean(IR_trial);
        stats{ipatient}.SI_trialavg{itemp}               = nanmean(SI_trial);
        stats{ipatient}.burstperc{itemp}                 = sum(stats{ipatient}.burst_trialsum{itemp}) / stats{ipatient}.spikecount{itemp} * 100;
        
    end
end

tbl = table;
iunit = 1;
for ipatient =  1 : 3
    for itemp = 1 : length(SpikeTrials{ipatient}.label)
        tbl.nodule(iunit)   = ipatient;
        tbl.unit(iunit)     = itemp;
        tbl.CV(iunit)       = stats{ipatient}.CV_trialavg{itemp};
        tbl.CV2(iunit)      = stats{ipatient}.CV2_trialavg{itemp};
        tbl.BI(iunit)       = stats{ipatient}.burstindex{itemp};
        tbl.BP(iunit)       = stats{ipatient}.burstperc{itemp};
        tbl.FR(iunit)       = stats{ipatient}.mean_freq{itemp};
        tbl.FF(iunit)       = stats{ipatient}.FF{itemp};
        tbl.template_pt(iunit) = SpikeRateStats{ipatient}.template_pt(itemp);
        tbl.template_tp(iunit) = SpikeRateStats{ipatient}.template_tp(itemp);
        iunit               = iunit + 1;
    end
end

filename = fullfile(config{ipatient}.datasavedir,'spikestats_revision.xlsx');
writetable(tbl,filename)
filename = fullfile(config{ipatient}.datasavedir,'spikestats_revision.csv');
writetable(tbl,filename)












%% extra spike stats cf Pierre
for ipatient =  1 : 3
    
    % continuous data
    for itemp = 1 : length(SpikeRateStats{ipatient}.isi)
        
        isi = SpikeRateStats{ipatient}.isi{itemp}/1000;
        SpikeRateStats{ipatient}.mean_freq(itemp)     = 1 / mean(isi);
        SpikeRateStats{ipatient}.median_freq(itemp)   = 1 / median(isi);
        SpikeRateStats{ipatient}.mean_isi(itemp)      = mean(isi);
        SpikeRateStats{ipatient}.stdev_isi(itemp)     = std(isi);
        SpikeRateStats{ipatient}.var_isi(itemp)       = std(isi)^2;
        SpikeRateStats{ipatient}.fano_isi(itemp)      = SpikeRateStats{ipatient}.var_isi(itemp)   / SpikeRateStats{ipatient}.mean_isi(itemp);
        SpikeRateStats{ipatient}.CV_isi(itemp)        = SpikeRateStats{ipatient}.stdev_isi(itemp) / SpikeRateStats{ipatient}.mean_isi(itemp);
        SpikeRateStats{ipatient}.burstindex(itemp)    = sum(SpikeRateStats{ipatient}.isi{itemp} < 5) / sum(SpikeRateStats{ipatient}.isi{itemp} < 100);
    end
    
    % per pattern
    for ipattern = 1 : size(SpikeRateStats{ipatient}.isi_pattern_all,2)
        
        for itemp = 1 : length(SpikeRateStats{ipatient}.isi_pattern_all{ipattern}.isi)
            
            trialavg_isi_all    = [];
            trialavg_isi_bl     = [];
            trialavg_isi_ac     = [];
            
            i = 1;
            for itrial                 = unique(SpikeRateStats{ipatient}.isi_pattern_all{ipattern}.trial{itemp})
                trialindx              = SpikeRateStats{ipatient}.isi_pattern_all{ipattern}.trial{itemp} == itrial;
                trialavg_isi_all(i)    = nanmean(SpikeRateStats{ipatient}.isi_pattern_all{ipattern}.isi{itemp}(trialindx));
                i = i + 1;
            end
            i = 1;
            for itrial                 = unique(SpikeRateStats{ipatient}.isi_pattern_bl{ipattern}.trial{itemp})
                trialindx              = SpikeRateStats{ipatient}.isi_pattern_bl{ipattern}.trial{itemp} == itrial;
                isi                    = SpikeRateStats{ipatient}.isi_pattern_bl{ipattern}.isi{itemp}(trialindx);     
                trialavg_isi_bl(i)     = nanmean(isi);               
                trail_CV_bl(i)         = nanstd(isi) / nanmean(isi);                
                x                      = isi(1:end-1) ./ isi(2:end);
                trial_CV2_bl(i)        = nanmean(2 * abs(x - 1) ./ (x + 1));
                trial_spikecount_bl(i) = sum(trialindx);                
                i = i + 1;
            end
            i = 1;
            trial_CV = [];
            for itrial                 = unique(SpikeRateStats{ipatient}.isi_pattern_ac{ipattern}.trial{itemp})
                trialindx              = SpikeRateStats{ipatient}.isi_pattern_ac{ipattern}.trial{itemp} == itrial;
                trialavg_isi_ac(i)     = nanmean(SpikeRateStats{ipatient}.isi_pattern_ac{ipattern}.isi{itemp}(trialindx));
                i = i + 1;
            end
            
            
            SpikeRateStats{ipatient}.freq_all(ipattern,itemp)   = 1 / nanmean(trialavg_isi_all);
            SpikeRateStats{ipatient}.freq_bl(ipattern,itemp)    = 1 / nanmean(trialavg_isi_bl);
            SpikeRateStats{ipatient}.freq_ac(ipattern,itemp)    = 1 / nanmean(trialavg_isi_ac);
            
            SpikeRateStats{ipatient}.freq_all(ipattern,itemp)   = 1 / nanmean(trialavg_isi_all);
            SpikeRateStats{ipatient}.freq_bl(ipattern,itemp)    = 1 / nanmean(trialavg_isi_bl);
            SpikeRateStats{ipatient}.freq_ac(ipattern,itemp)    = 1 / nanmean(trialavg_isi_ac);
            
%             SpikeRateStats{ipatient}.freq_bl(ipattern,itemp)    = stats_bar{ipatient}.clusterstat{ilabel}{itemp}.bl.avg;  
            
            SpikeRateStats{ipatient}.fano_all(ipattern,itemp)   = nanstd(trialavg_isi_all)^2 / nanmean(trialavg_isi_all);
% %             SpikeRateStats{ipatient}.fano_bl(ipattern,itemp)    = nanstd(trialavg_isi_bl)^2  / nanmean(trialavg_isi_bl);
            SpikeRateStats{ipatient}.fano_ac(ipattern,itemp)    = nanstd(trialavg_isi_ac)^2  / nanmean(trialavg_isi_ac);
            SpikeRateStats{ipatient}.cv_all(ipattern,itemp)     = nanstd(trialavg_isi_all)   / nanmean(trialavg_isi_all);
%             SpikeRateStats{ipatient}.cv_bl(ipattern,itemp)      = nanstd(trialavg_isi_bl)    / nanmean(trialavg_isi_bl);
            SpikeRateStats{ipatient}.cv_ac(ipattern,itemp)      = nanstd(trialavg_isi_ac)    / nanmean(trialavg_isi_ac);
            SpikeRateStats{ipatient}.BI_all(ipattern,itemp)     = sum(SpikeRateStats{ipatient}.isi_pattern_all{ipattern}.isi{itemp} < 0.005) / sum(SpikeRateStats{ipatient}.isi_pattern_all{ipattern}.isi{itemp} < 0.100);
            SpikeRateStats{ipatient}.BI_bl(ipattern,itemp)      = sum(SpikeRateStats{ipatient}.isi_pattern_bl{ipattern}.isi{itemp} < 0.005)  / sum(SpikeRateStats{ipatient}.isi_pattern_bl{ipattern}.isi{itemp} < 0.100);
            SpikeRateStats{ipatient}.BI_ac(ipattern,itemp)      = sum(SpikeRateStats{ipatient}.isi_pattern_ac{ipattern}.isi{itemp} < 0.005)  / sum(SpikeRateStats{ipatient}.isi_pattern_ac{ipattern}.isi{itemp} < 0.100);
            
            SpikeRateStats{ipatient}.cv_bl(ipattern,itemp)      = nanmean(trail_CV_bl);
            SpikeRateStats{ipatient}.cv2_bl(ipattern,itemp)     = nanmean(trial_CV2_bl);
            SpikeRateStats{ipatient}.fano_bl(ipattern,itemp)    = nanstd(trial_spikecount_bl)^2  / nanmean(trial_spikecount_bl);     
        end
    end
end


tbl = table;
iunit = 1;
for ipatient =  1 : 3
    
    for itemp = 1 : length(SpikeRateStats{ipatient}.isi_pattern_all{ipattern}.isi)
        
        
        tbl.nodule(iunit)   = ipatient;
        total               = length(SpikeRateStats{ipatient}.isi{itemp});
        short               = sum(SpikeRateStats{ipatient}.isi{itemp} < 2);
        long                = sum(SpikeRateStats{ipatient}.isi{itemp} >= 2);
        tbl.unit(iunit)     = itemp;
        tbl.total(iunit)    = total;
        
        if short/total * 100 < 1
            tbl.percRPV{iunit} = round(short/total * 100,2);
        else
            tbl.percRPV{iunit} = round(short/total * 100,2);
        end
        
        tbl.template_pt(iunit)      = SpikeRateStats{ipatient}.template_pt(itemp);
        tbl.template_tp(iunit)      = SpikeRateStats{ipatient}.template_tp(itemp);        
        tbl.template_width(iunit)   = SpikeRateStats{ipatient}.template_width(itemp);
   
        tbl.ipatient(iunit)         = ipatient;
        tbl.ALL_mean_freq(iunit)    = SpikeRateStats{ipatient}.mean_freq(itemp);
        tbl.ALL_median_freq(iunit)  = SpikeRateStats{ipatient}.median_freq(itemp);
        tbl.ALL_fano(iunit)         = SpikeRateStats{ipatient}.fano_isi(itemp);
        tbl.ALL_CV(iunit)           = SpikeRateStats{ipatient}.CV_isi(itemp);
        tbl.ALL_BI(iunit)           = SpikeRateStats{ipatient}.burstindex(itemp);
        
        if ipatient ~= 3
            
            tbl.SW_bl_freq(iunit)   = SpikeRateStats{ipatient}.freq_bl(1,itemp);
            tbl.FA_bl_freq(iunit)   = SpikeRateStats{ipatient}.freq_bl(2,itemp);
            tbl.ES_bl_freq(iunit)   = SpikeRateStats{ipatient}.freq_bl(3,itemp);
            tbl.SW_ac_freq(iunit)   = SpikeRateStats{ipatient}.freq_ac(1,itemp);
            tbl.FA_ac_freq(iunit)   = SpikeRateStats{ipatient}.freq_ac(2,itemp);
            tbl.ES_ac_freq(iunit)   = SpikeRateStats{ipatient}.freq_ac(3,itemp);
            
            tbl.SW_bl_fano(iunit)   = SpikeRateStats{ipatient}.fano_bl(1,itemp);
            tbl.FA_bl_fano(iunit)   = SpikeRateStats{ipatient}.fano_bl(2,itemp);
            tbl.ES_bl_fano(iunit)   = SpikeRateStats{ipatient}.fano_bl(3,itemp);
            tbl.SW_ac_fano(iunit)   = SpikeRateStats{ipatient}.fano_ac(1,itemp);
            tbl.FA_ac_fano(iunit)   = SpikeRateStats{ipatient}.fano_ac(2,itemp);
            tbl.ES_ac_fano(iunit)   = SpikeRateStats{ipatient}.fano_ac(3,itemp);
            
            tbl.SW_bl_CV(iunit)     = SpikeRateStats{ipatient}.cv_bl(1,itemp);
            tbl.FA_bl_CV(iunit)     = SpikeRateStats{ipatient}.cv_bl(2,itemp);
            tbl.ES_bl_CV(iunit)     = SpikeRateStats{ipatient}.cv_bl(3,itemp);
            tbl.SW_ac_CV(iunit)     = SpikeRateStats{ipatient}.cv_ac(1,itemp);
            tbl.FA_ac_CV(iunit)     = SpikeRateStats{ipatient}.cv_ac(2,itemp);
            tbl.ES_ac_CV(iunit)     = SpikeRateStats{ipatient}.cv_ac(3,itemp);
            
            tbl.SW_bl_CV2(iunit)    = SpikeRateStats{ipatient}.cv2_bl(1,itemp);
            tbl.FA_bl_CV2(iunit)    = SpikeRateStats{ipatient}.cv2_bl(2,itemp);
            tbl.ES_bl_CV2(iunit)    = SpikeRateStats{ipatient}.cv2_bl(3,itemp);
            
            tbl.SW_bl_BI(iunit)     = SpikeRateStats{ipatient}.BI_bl(1,itemp);
            tbl.FA_bl_BI(iunit)     = SpikeRateStats{ipatient}.BI_bl(2,itemp);
            tbl.ES_bl_BI(iunit)     = SpikeRateStats{ipatient}.BI_bl(3,itemp);
            tbl.SW_ac_BI(iunit)     = SpikeRateStats{ipatient}.BI_ac(1,itemp);
            tbl.FA_ac_BI(iunit)     = SpikeRateStats{ipatient}.BI_ac(2,itemp);
            tbl.ES_ac_BI(iunit)     = SpikeRateStats{ipatient}.BI_ac(3,itemp);
            
        else
            
            tbl.FA_bl_freq(iunit)   = SpikeRateStats{ipatient}.freq_bl(1,itemp);
            tbl.ES_bl_freq(iunit)   = SpikeRateStats{ipatient}.freq_bl(2,itemp);
            tbl.FA_ac_freq(iunit)   = SpikeRateStats{ipatient}.freq_ac(1,itemp);
            tbl.ES_ac_freq(iunit)   = SpikeRateStats{ipatient}.freq_ac(2,itemp);
            
            tbl.FA_bl_fano(iunit)   = SpikeRateStats{ipatient}.fano_bl(1,itemp);
            tbl.ES_bl_fano(iunit)   = SpikeRateStats{ipatient}.fano_bl(2,itemp);
            tbl.FA_ac_fano(iunit)   = SpikeRateStats{ipatient}.fano_ac(1,itemp);
            tbl.ES_ac_fano(iunit)   = SpikeRateStats{ipatient}.fano_ac(2,itemp);
            
            tbl.FA_bl_CV(iunit)     = SpikeRateStats{ipatient}.cv_bl(1,itemp);
            tbl.ES_bl_CV(iunit)     = SpikeRateStats{ipatient}.cv_bl(2,itemp);
            
            tbl.FA_bl_CV2(iunit)    = SpikeRateStats{ipatient}.cv2_bl(1,itemp);
            tbl.ES_bl_CV2(iunit)    = SpikeRateStats{ipatient}.cv2_bl(2,itemp);
            
            tbl.FA_ac_CV(iunit)     = SpikeRateStats{ipatient}.cv_ac(1,itemp);
            tbl.ES_ac_CV(iunit)     = SpikeRateStats{ipatient}.cv_ac(2,itemp);
            
            tbl.FA_bl_BI(iunit)     = SpikeRateStats{ipatient}.BI_bl(1,itemp);
            tbl.ES_bl_BI(iunit)     = SpikeRateStats{ipatient}.BI_bl(2,itemp);
            tbl.FA_ac_BI(iunit)     = SpikeRateStats{ipatient}.BI_ac(1,itemp);
            tbl.ES_ac_BI(iunit)     = SpikeRateStats{ipatient}.BI_ac(2,itemp);
        end
        
        iunit                       = iunit + 1;
        
    end
end


filename = fullfile(config{ipatient}.datasavedir,'spikestats_revision.xlsx');
writetable(tbl,filename)
filename = fullfile(config{ipatient}.datasavedir,'spikestats_revision.csv');
writetable(tbl,filename)






%     temp                = dir(fullfile(config{ipatient}.datasavedir,[config{ipatient}.prefix,'all_data_',config{ipatient}.circus.channel{1}(1:end-2),'_*.ncs']));
%     hdr_fname           = fullfile(temp(1).folder,temp(1).name);
%     hdr                 = ft_read_header(hdr_fname); % take the first file to extract the header of the data
%     
%     % redefine trials to 1-second windows for ISI
%     cfgtemp             = [];
%     cfgtemp.trl         = (1 : hdr.Fs*1 : hdr.nSamples)';
%     cfgtemp.trl(:,2)    = cfgtemp.trl(:,1) + hdr.Fs*1 - 1;
%     cfgtemp.trl(:,3)    = zeros(size(cfgtemp.trl,1),1);
%     cfgtemp.trl         = cfgtemp.trl(1:end-1,:);
%     cfgtemp.trlunit     = 'samples';
%     cfgtemp.hdr         = hdr;
%     spiketrials_1s      = ft_spike_maketrials(cfgtemp,SpikeRaw);
%     
%     for itemp = 1 : length(SpikeRateStats.isi)
%         
%         for itrial = unique(spiketrials_1s.trial{itemp})
%             trialsel = find(spiketrials_1s.trial{itemp} == itrial);
%             if trialsel > 1
%                 temp(itemp,itrial) = mean( diff( spiketrials_1s.time{itemp}(trialsel)));
%             else
%                 temp(itemp,itrial) = nan;
%             end
%         end
%         
%     end
%     temp(temp == 0) = nan;
%     t = nanmean(temp,2)';
% 
% end

%% average spike LFP averages

w_pre  = config{ipatient}.spike.pre * spike_LFP{itemp}.fsample;
w_post = config{ipatient}.spike.post * spike_LFP{itemp}.fsample;

for itemp = 1 : size(spike_LFP,2)
    
    spike_LFP{itemp}.trial_norm = [];
    
    % average
    i = 1;
    for trialnr = 1 : size(spike_LFP{itemp}.trial,2)
        if size(find(SpikeRaw.sample{itemp} < SpikeRaw.sample{itemp}(trialnr)+w_post & SpikeRaw.sample{itemp} > SpikeRaw.sample{itemp}(trialnr)-w_pre),1) == 1
            
            fprintf('Normalizing unit %d of %d: trial %d of %d\n',itemp,size(spike_LFP,2),trialnr,size(spike_LFP{itemp}.trial,2));
            
            %             spike_LFP{itemp}.trial_norm(i,:) = (spike_LFP{itemp}.trial{trialnr} - min(spike_LFP{itemp}.trial{trialnr})) / (max(spike_LFP{itemp}.trial{trialnr}) - min(spike_LFP{itemp}.trial{trialnr}) );
            spike_LFP{itemp}.trial_norm{i} = spike_LFP{itemp}.trial{trialnr};
            
            i = i + 1;
        else
            fprintf('Skipping unit %d of %d: trial %d of %d because it overlaps with other spikes\n',itemp,size(spike_LFP,2),trialnr,size(spike_LFP{itemp}.trial,2));
        end
        
    end
end

for itemp = 1 : size(spike_LFP,2)
    cfgtemp = [];
    cfgtemp.demean = 'yes';
    cfgtemp.baselinewindow = [-0.001, -0.0005];
    spike_LFP{itemp} = ft_preprocessing(cfgtemp,spike_LFP{itemp});
    spike_LFP_avg{itemp} = ft_timelockanalysis([],spike_LFP{itemp});
    spike_LFP_avg{itemp}.avg_norm = (spike_LFP_avg{itemp}.avg - min(spike_LFP_avg{itemp}.avg)) / (max(spike_LFP_avg{itemp}.avg) - min(spike_LFP_avg{itemp}.avg) );
end


figure; hold;
for itemp = 1 : size(spike_LFP,2)
    plot(spike_LFP_avg{itemp}.time,spike_LFP_avg{itemp}.avg_norm,'k');
end


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
        xi = xi + 1;
        
        for icontact = 1 : length(stat_TFR_micro{ipatient}{imarker}.corrs.avg)
            xpos = (xi-1) * 3;;
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

