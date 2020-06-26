function plotTFRcontinuous(cfg, TFR)

  % PLOTTFRCONTINUOUS creates statistics based on hypnogram and markers in MuseStruct
  %
  % use as
  %   plotTFRcontinuous(cfg, TFR)

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

for ipart = 1 : size(TFR,2)
    fig = figure;

    for ichan = 1 : size(TFR{ipart}.label,1)

        subplot(size(TFR{ipart}.label,1),1,ichan);

        cfgtemp                 = [];
        cfgtemp.channel         = TFR{ipart}.label{ichan};
        %         cfgtemp.ylim            = [1 40];
        cfgtemp.baseline        = [TFR{1}.time(1)+120 TFR{1}.time(end)-120];
        cfgtemp.baselinetype    = 'relative';
        cfgtemp.colorbar        = 'no';
        cfgtemp.colorbar        = 'yes';
        %     cfgtemp.zlim            = 'maxabs';
        %     cfgtemp.xlim            = config{ipatient}.;
        %     cfgtemp.title           = 'Relative change from Baseline';
        cfgtemp.parameter       = 'powspctrm';
        cfgtemp.colormap        = parula(5000);
        cfgtemp.renderer        = 'painters';
        ft_singleplotTFR(cfgtemp,TFR{ipart});
    end

    % print to file
    set(fig,'PaperOrientation','landscape');
    set(fig,'PaperUnits','normalized');
    set(fig,'PaperPosition', [0 0 1 1]);
    print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'p',num2str(ipart),'-TFR_',TFR{ipart}.label{ichan},'.pdf']),'-r600');
    close all
end
