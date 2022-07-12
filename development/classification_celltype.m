function celltype = classification_celltype(config, force, do_plot)
% 
% Use as 
%       celltype = classification_celltype(config, force, do_plot)
% 
% force   : whether to redo analyses or read previous save (true/false)
% do_plot : optional (default = true), can be set to false to only compute
%           value without plotting
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

%% The 2 references used to make this function:
% 
% Peyrache et al., 2012 - Spatiotemporal dynamics of neocortical excitation and
% inhibition during human sleep. www.pnas.org/cgi/doi/10.1073/pnas.1109895109
% Human recordings
% Discrimination of Pyramidal (Pyr) Cells and Fast-Spiking (FS) Interneurons (Int).
% Average waveforms were computed for each isolated cell. As described previously,
% the half width of the extracellular positive deflection has, at the neuronal population level,
% a bimodal distribution (1, 2). The separation is even more
% striking when the valley-to-peak parameter (2) is added for 2D
% clustering (Fig. 2 A and B). Automatic clustering of these average waveforms
% from individual cells by using a k-means algorithm discriminated two groups
% of cells (Fig. 2 A and B). The resultant clustering was further confirmed by
% an E-M clustering method (Fig. S1)

% Elahian et al., 2018 - Low-Voltage Fast Seizures in Humans Begin with Increased
% Interneuron Firing. doi:10.1002/ana.25325
% Utilise le peak asymetry en plus, defini comme :
% peak_asymetry = (b-a)/(b+a);
% With a : amplitude of the first peak
% And b : amplitude of the second peak
% (with the AP trough between the two peaks)


%% analysis

config{1}.tablesavedir = ft_getopt(config{1}, 'tablesavedir', config{1}.datasavedir);
fname = fullfile(config{1}.tablesavedir, 'allpatients_cell_type.xlsx');

if exist(fname, 'file') && force == false
    fprintf('Reading %s\n', fname);
    celltype = readtable(fname);
    return
end

if nargin < 3
    do_plot = true;
end

%gather all data to perform clustering
data.patient_nr = [];
data.ipart      = [];
data.label      = {};
data.cluster_group = {};
data.halfwidth  = [];
data.peaktrough = [];
data.troughpeak = [];
data.peak_asymetry = [];

for ipatient = 1:size(config, 2)
    waveformstats = spikeWaveformStats(config{ipatient}, [], false);
    for ipart = 1:size(config{ipatient}.directorylist, 2)
        for i_unit = 1:size(waveformstats{ipart}.label, 2)
            data.patient_nr(end+1) = ipatient;
            data.ipart(end+1)      = ipart;
            data.label{end+1}      = waveformstats{ipart}.label{i_unit};
            data.cluster_group{end+1} = strrep(waveformstats{ipart}.cluster_group{i_unit}, ' ', '');
            data.halfwidth(end+1)  = waveformstats{ipart}.halfwidth.val(i_unit)*1000;%ms
            data.troughpeak(end+1) = waveformstats{ipart}.troughpeak.val(i_unit) * 1000;
            
            a = waveformstats{ipart}.peaktrough.y(i_unit, 1) * -waveformstats{ipart}.peak_direction(i_unit);
            b = waveformstats{ipart}.troughpeak.y(i_unit, 2) * -waveformstats{ipart}.peak_direction(i_unit);
            data.peak_asymetry(end+1) = (b-a)/(b+a);
        end
    end
end

%remove nans:
toremove = isnan(data.halfwidth) | isnan(data.troughpeak);
data.halfwidth(toremove)     = [];
data.troughpeak(toremove)    = [];
data.peak_asymetry(toremove) = [];
data.label(toremove)         = [];
data.cluster_group(toremove) = [];

%clusters peyrache
X                      = [data.halfwidth; data.troughpeak]';
[idx_p,C_p]            = kmeans(X, 2, 'maxiter', 100000, 'replicates', 100);
[~, idx_in]            = min([mean(data.halfwidth(idx_p == 1)), mean(data.halfwidth(idx_p == 2))]);
idx_celltype.p         = ["putative_PN", "putative_PN"];
idx_celltype.p(idx_in) = "putative_IN";

%clusters elahian
X                      = [data.halfwidth; data.troughpeak; data.peak_asymetry]';
[idx_e,C_e]            = kmeans(X, 2, 'maxiter', 100000, 'replicates', 100);
[~, idx_in]            = min([mean(data.halfwidth(idx_e == 1)), mean(data.halfwidth(idx_e == 2))]);
idx_celltype.e         = ["putative_PN", "putative_PN"];
idx_celltype.e(idx_in) = "putative_IN";

% table cell types
celltype = table.empty;
irow = 0;
for i_unit = 1:size(data.label, 2)
    irow = irow+1;
    celltype.patient_ID{irow}         = config{data.patient_nr(i_unit)}.prefix(1:end-1);
    celltype.unit_ID{irow}            = data.label{i_unit};
    celltype.cluster_group{irow}      = data.cluster_group{i_unit};
    celltype.halfwidth_ms(irow)       = data.halfwidth(i_unit);
    celltype.troughpeak_ms(irow)      = data.troughpeak(i_unit);
    celltype.peak_asymetry(irow)      = data.peak_asymetry(i_unit);
    celltype.celltype_peyrache{irow}  = idx_celltype.p(idx_p(i_unit));
    celltype.celltype_elahian{irow}   = idx_celltype.e(idx_e(i_unit));
end

%% plots
if do_plot
    for imethod = ["peyrache", "elahian"]
        fig = figure;
        
        for i_unit = 1:size(data.label, 2)
            
            celltype_unit = celltype.(sprintf('celltype_%s', imethod)){i_unit};
            if contains(celltype_unit, 'PN')
                plottype = 'b';
            elseif contains(celltype_unit, 'IN')
                plottype = 'r';
            end
            
            x = data.halfwidth(i_unit);
            y = data.troughpeak(i_unit); %en ms
            peak_asymetry = data.peak_asymetry(i_unit);
            
            if imethod == "peyrache"
                if contains(data.cluster_group{i_unit}, 'good')
                    scatter(x,y,plottype, 'filled');
                else
                    scatter(x,y,plottype);
                end
                %plot centroids
                if i_unit == size(data.label, 2)
                    scatter(C_p(:,1), C_p(:,2), 'xk', 'linewidth', 2);
                end
            elseif imethod == "elahian"
                if contains(data.cluster_group{i_unit}, 'good')
                    scatter3(x,y,peak_asymetry,plottype, 'filled'); % 3D plot
                else
                    scatter3(x,y,peak_asymetry,plottype); % 3D plot
                end
                %plot centroids
                if i_unit == size(data.label, 2)
                    scatter3(C_e(:,1), C_e(:,2), C_e(:,3), 'xk', 'linewidth', 2);
                end
            end
            hold on;
        end
        
        
        set(gca, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 15);
        xlabel('halfwidth(ms)');
        ylabel('trough-peak (ms)');
        zlabel('peak asymetry');
        
        figname = fullfile(config{1}.imagesavedir, 'Classification_in_pn', imethod);
        savefigure_own(fig, figname, 'png', 'pdf', 'fig', 'close');
        
    end
end

%% save results
delete(fname);
writetable(celltype, fname);

