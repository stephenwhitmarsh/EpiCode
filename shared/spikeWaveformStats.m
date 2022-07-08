function 	stats = spikeWaveformStats(cfg, SpikeWaveforms, force)

% SPIKEWAVEFORMSTATS computes average and variance of the spike waveform.
% Then, the average is interpolated and the following parameters are
% computed : amplitude, halfwidth, peaktrough and trouhpeak.
% Input spike waveform data is get from readSpikeWaveforms.m.
%
% Use as :
%   stats = spikeWaveformStats(cfg, SpikeWaveforms, force)
%
% Note :
% - spikes must be aligned so they can be averaged
% - spikes can be positive or negative
%
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


fname_out = fullfile(cfg.datasavedir, strcat(cfg.prefix, 'spikewaveformstats.mat'));

if exist(fname_out, 'file') && force == false
    fprintf('Loading precomputed spike waveform stats\n');
    load(fname_out, 'stats');
    return
else
    fprintf('(re-) computing spike waveform stats\n');
end

% add external/intersections path if not already
mfile_name = mfilename('fullpath');
pathstr    = fileparts(fileparts(mfile_name));
pathCell   = regexp(path, pathsep, 'split');
if ispc  % Windows is not case-sensitive
    onPath = any(strcmpi(fullfile(pathstr, ['external', filesep, 'intersections']), pathCell));
else
    onPath = any(strcmp(fullfile(pathstr, ['external', filesep, 'intersections']), pathCell));
end
if ~onPath
    addpath(fullfile(pathstr, ['external', filesep, 'intersections']));
end

for ipart = 1:size(SpikeWaveforms)
        ft_progress('init', 'text');
        for icluster = 1:size(SpikeWaveforms{ipart}, 2)

            ft_progress(0, 'Compute waveform for part %d, cluster %d/%d', ipart, icluster, size(SpikeWaveforms{ipart}, 2));
            ok = true;

            if ~isempty(SpikeWaveforms{ipart}{icluster})

                %average waveform
                evalc('waveformavg = ft_timelockanalysis([], SpikeWaveforms{ipart}{icluster});');%remove text output

                %interpolate the averaged waveform to have more precise values
                time_temp        = linspace(waveformavg.time(1),waveformavg.time(end),1000);
                avg_temp         = pchip(waveformavg.time,waveformavg.avg,time_temp);
                var_temp         = pchip(waveformavg.time,waveformavg.var,time_temp);
                waveformavg.time = time_temp;
                waveformavg.avg  = avg_temp;
                waveformavg.var  = var_temp;
                %scatter(waveformavg.time, waveformavg.avg, '.');

                %search if AP is positive or negative
                flip = sign(min(waveformavg.avg) + max(waveformavg.avg));

                %amplitude
                bl               = min(waveformavg.avg.*flip);
                [peak, peak_idx] = max(waveformavg.avg.*flip);
                amplitude.val    = peak-bl;
                amplitude.x      = waveformavg.time(peak_idx);
                amplitude.y      = waveformavg.avg(peak_idx);
                %scatter(amplitude.x, amplitude.y, 'xk');

                %halfwidth
                halfamp     = bl+(peak-bl)/2;
                x  = waveformavg.time;
                y1 = waveformavg.avg;
                y2 = ones(size(x)) .* halfamp .* flip;
                [x_intersect, y_intersect] = intersections(x,y1,x,y2,true);
                
                if length(x_intersect) >= 2 && any(x_intersect < amplitude.x) && any(x_intersect > amplitude.x)
                  
                    idx = find(x_intersect < amplitude.x, 1, 'last');
                    halfwidth.val = diff(x_intersect([idx, idx+1]));
                    halfwidth.x   = x_intersect([idx, idx+1]);
                    halfwidth.y   = y_intersect([idx, idx+1]);
                    %scatter(halfwidth.x, halfwidth.y, 'xk');

                    % Find throughs
                    [Yneg,Xneg_temp] = findpeaks(-waveformavg.avg.*flip,waveformavg.time);
                    
                    if length(Xneg_temp) >= 2 && any(Xneg_temp < amplitude.x) && any(Xneg_temp > amplitude.x)
                        % Search first through before and after peak
                        [Xneg(1),x_idx] = max(Xneg_temp(Xneg_temp < amplitude.x));
                        Xneg(2)         = min(Xneg_temp(Xneg_temp > amplitude.x));
                        Yneg            = Yneg([x_idx, x_idx+1]);

                        peaktrough.val  = abs(amplitude.x-Xneg(1));
                        peaktrough.x    = [Xneg(1)  amplitude.x];
                        peaktrough.y    = [-Yneg(1)*flip amplitude.y];
                        troughpeak.val  = abs(amplitude.x-Xneg(2));
                        troughpeak.x    = [amplitude.x Xneg(2)];
                        troughpeak.y    = [amplitude.y -Yneg(2)*flip];
                        %scatter(peaktrough.x, peaktrough.y, 'xk');
                        %scatter(troughpeak.x, troughpeak.y, 'xk');
                        
                    else
                        ok = false;
                    end
                else
                    ok = false;
                end
            else
                ok = false;
            end %isempty

            %store for output
            
            if isempty(SpikeWaveforms{ipart}{icluster})
                stats{ipart}.label{icluster}          = [];
                stats{ipart}.waveformavg{icluster}    = [];
                stats{ipart}.peak_direction(icluster) = nan;
            else
                stats{ipart}.label{icluster}          = SpikeWaveforms{ipart}{icluster}.label{1};
                stats{ipart}.cluster_group{icluster}  = SpikeWaveforms{ipart}{icluster}.cluster_group;
                stats{ipart}.waveformavg{icluster}    = waveformavg;
                stats{ipart}.peak_direction(icluster) = flip;
            end
            if ok
                stats{ipart}.amplitude.val(icluster)  = amplitude.val;
                stats{ipart}.amplitude.x(icluster)    = amplitude.x;
                stats{ipart}.amplitude.y(icluster)    = amplitude.y;
                stats{ipart}.halfwidth.val(icluster)  = halfwidth.val;
                stats{ipart}.halfwidth.x(icluster,:)  = halfwidth.x;
                stats{ipart}.halfwidth.y(icluster,:)  = halfwidth.y;
                stats{ipart}.peaktrough.val(icluster) = peaktrough.val;
                stats{ipart}.peaktrough.x(icluster,:) = peaktrough.x;
                stats{ipart}.peaktrough.y(icluster,:) = peaktrough.y;
                stats{ipart}.troughpeak.val(icluster) = troughpeak.val;
                stats{ipart}.troughpeak.x(icluster,:) = troughpeak.x;
                stats{ipart}.troughpeak.y(icluster,:) = troughpeak.y;
            else
                stats{ipart}.amplitude.val(icluster)  = nan;
                stats{ipart}.amplitude.x(icluster)    = nan;
                stats{ipart}.amplitude.y(icluster)    = nan;
                stats{ipart}.halfwidth.val(icluster)  = nan;
                stats{ipart}.halfwidth.x(icluster,:)  = [nan nan];
                stats{ipart}.halfwidth.y(icluster,:)  = [nan nan];
                stats{ipart}.peaktrough.val(icluster) = nan;
                stats{ipart}.peaktrough.x(icluster,:) = [nan nan];
                stats{ipart}.peaktrough.y(icluster,:) = [nan nan];
                stats{ipart}.troughpeak.val(icluster) = nan;
                stats{ipart}.troughpeak.x(icluster,:) = [nan nan];
                stats{ipart}.troughpeak.y(icluster,:) = [nan nan];
            end
        end
        ft_progress('close');
%     end
end

save(fname_out, 'stats', '-v7.3');
