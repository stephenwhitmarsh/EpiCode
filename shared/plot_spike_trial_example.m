function plot_spike_trial_example(cfg, spiketrials, waveformstats, ipart, markername, itrial, maxsize)

% Use as : 
%       plot_spike_trial_example(cfg, spiketrials, waveformstats, ipart, markername, itrial, maxsize)
% plot_spike_trial_example plots, for one trial, the raw LFP, the LFP
% high-pass filtered, and the LFP high-pass filtered colored at the moments
% where Spiking-Circus found spikes.

% Input: 
% - spiketrials   : output from readSpikeTrials.m
% - waveformstats : output from readSpikeWaveforms.m
% - ipart         : nr. of the part to use
% - markername    : name of the analysis to use
% - itrial        : number of the trial to use 
% - maxsize       : (optional, default = inf) if trial is longer than maxsize, then 
%                    the first maxsize seconds are displayed

% This file is part of EpiCode, see
% http://www.github.com/stephenwhitmarsh/EpiCode for documentation and details.
%
%   EpiCode is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   EpiCode is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with EpiCode. If not, see <http://www.gnu.org/licenses/>.
% 

if nargin < 7
    maxsize = inf;
end

for ichannel = 1:size(cfg.circus.channel, 2)
    
    %find trial samples to load filtered data
    idir     = spiketrials{ipart}.(markername).trialinfo.idir(itrial);
    channame = cfg.circus.channel{ichannel};
    temp     = dir(fullfile(cfg.rawdir, cfg.directorylist{ipart}{idir}, sprintf('*%s.ncs', channame)));
    if strcmp(channame(1), '_')
        channame = channame(2:end);
    end
    datapath = fullfile(temp.folder, temp.name);
    Fs       = spiketrials{ipart}.(markername).hdr.Fs;
    
    cfgtemp             = [];
    cfgtemp.trl(1)      = spiketrials{ipart}.(markername).trialinfo.begsample_dir(itrial);
    cfgtemp.trl(2)      = cfgtemp.trl(1) + (spiketrials{ipart}.(markername).trialinfo.endsample(itrial) - spiketrials{ipart}.(markername).trialinfo.begsample(itrial));
    cfgtemp.trl(3)      = spiketrials{ipart}.(markername).trialinfo.offset(itrial);
    cfgtemp.dataset     = datapath;
    
    %shorten trial if it is too long
    if cfgtemp.trl(2) - cfgtemp.trl(1) > maxsize * Fs
        cfgtemp.trl(2) = cfgtemp.trl(1) + maxsize * Fs;
    end
    cfgtemp.bsfilter    = 'yes';
    cfgtemp.bsfreq      = [49 51];
    cfgtemp.bsinstabilityfix = 'reduce';
    dat                 = ft_preprocessing(cfgtemp);
    
    cfgtemp.hpfilter    = 'yes';
    cfgtemp.hpfreq      = 300;
    cfgtemp.hpfiltord   = 3;
    dat_hpfilt          = ft_preprocessing(cfgtemp);
    
    iplot = 1;
    fig = figure;
    sgtitle(sprintf('%s %s (%s, trial %d)', cfg.prefix(1:end-1), markername, channame, itrial),...
        'interpreter', 'none', 'FontWeight', 'bold', 'Fontsize', 18);
    
    %plot LFP
    subplot(3,1,iplot); hold on;
    flip = -waveformstats{ipart}.peak_direction(1);
    p = plot(dat.time{1}, dat.trial{1}*flip, 'k', 'linewidth', 1);
    l1 = legend(p, channame, 'location', 'eastoutside', 'interpreter', 'none');
    set(gca, 'tickdir', 'out', 'fontsize', 15);
    ylabel('uV');
    axis tight
    
    %plot MUA without and with higlights
    for do_highlight = [false true]
        
        iplot = iplot +1;
        %plot data
        subplot(3,1,iplot); hold on;
        flip = -waveformstats{ipart}.peak_direction(1);
        p = plot(dat_hpfilt.time{1}, dat_hpfilt.trial{1}*flip, 'k', 'linewidth', 1);
        
        if do_highlight
            unitssel = ichannel == spiketrials{ipart}.(markername).template_maxchan +1;
            C = linspecer(sum(unitssel));
            icolor = 0;
            
            ileg = 0;
            leg = [];
            for i_unit = find(unitssel)
                icolor = icolor+1;
                idx = spiketrials{ipart}.(markername).trial{i_unit} == itrial;
                spikes = spiketrials{ipart}.(markername).time{i_unit}(idx);
                first_spike = true;
                for ispike = 1 : size(spikes, 2)
                    if spikes(ispike) - 0.02 < dat_hpfilt.time{1}(1) || spikes(ispike) + 0.02 > dat_hpfilt.time{1}(end)
                        continue
                    end
                    if first_spike
                        ileg = ileg +1;
                        first_spike = false;
                    end
                    temp = round((spikes(ispike) - 0.001) * Fs) : round((spikes(ispike) + 0.02) * Fs);
                    sel = double(temp) - double(spiketrials{ipart}.(markername).trialinfo.offset(itrial));
                    leg{ileg} = plot(dat_hpfilt.time{1}(sel), dat_hpfilt.trial{1}(sel)*flip, 'color', C(icolor, :), 'linewidth', 1.5);%[1 0.6 0.2]);
                    leg_name{ileg} = spiketrials{ipart}.(markername).label{i_unit};
                end
            end
            if ~isempty(leg)
                l3 = legend([leg{:}], leg_name{:}, 'Location', 'EastOutside', 'interpreter', 'none');
            else
                l3 = legend(p, channame, 'Location', 'EastOutside', 'interpreter', 'none');
            end
        else
            l2 = legend(p, channame, 'Location', 'EastOutside', 'interpreter', 'none');
        end
        
        set(gca, 'tickdir', 'out', 'fontsize', 15);
        ylabel('uV');
        axis tight
        
    end
    
    xlabel('time (s)');
    
    %make legend box the same size to have aligned x axis
    maxlen = max(cellfun(@length, l3.String));
    difflen = maxlen - length(l1.String{1});
    l1.String{1}   = [l1.String{1}, repmat(' ', 1, difflen+3), '.'];
    l2.String{1}   = l1.String{1};
    l1.FontSize    = l3.FontSize;
    l2.FontSize    = l3.FontSize;
    
    figname = fullfile(cfg.imagesavedir, 'example_trial_spikes', sprintf('%s%s_%s_trial%d', cfg.prefix, channame, markername, itrial));
    savefigure_own(fig, figname, 'landscape', 'png', 'pdf', 'close');
    
end