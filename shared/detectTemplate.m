function [MuseStruct, Tindx_unique, LFP_avg] = detectTemplate(cfg, MuseStruct, template, force)

% DETECTEMPLATE detect templates using normalized crosscorrelation.
% Compares templates, and writes markers of most fitting template to Muse
% event file.
%
% use as
%    [MuseStruct, C_norm, Tindx_unique, LFP_avg] = detectTemplate(cfg, MuseStruct, template, force)
%
% Necessary fields (as defined in _setparams function):
%
% This file is part of EpiCode, see
% http://www.github.com/stephenwhitmarsh/EpiCode for documentation and details.
%
%    EpiCode is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    EpiCode is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with EpiCode. If not, see <http://www.gnu.org/licenses/>.

cfg.visible                 = ft_getopt(cfg, 'visible', 'on');
cfg.template.reref          = ft_getopt(cfg.template, 'reref', 'no');
cfg.template.refmethod      = ft_getopt(cfg.template, 'refmethod', 'none');
cfg.template.latency        = ft_getopt(cfg.template, 'latency', 'all');
cfg.template.writemuse      = ft_getopt(cfg.template, 'writemuse', false);
cfg.template.threshold      = ft_getopt(cfg.template, 'threshold', 4);
cfg.template.name           = ft_getopt(cfg.template, 'name', 'TemplateDetect');

fname_out                   = fullfile(cfg.datasavedir, [cfg.prefix, 'MuseStruct_detectedTemplates.mat']);

if nargin == 1
    if exist(fname_out, 'file')
        fprintf('Loading results template detection: %s\n', fname_out);
        % repeat to deal with load errors
        count = 0;
        err_count = 0;
        while count == err_count
            try
                load(fname_out, 'MuseStruct', 'Tindx_unique', 'LFP_avg');
            catch ME
                err_count = err_count + 1;
                disp('Something went wrong loading the file. Trying again...')
            end
            count = count + 1;
        end
        return;
    else
        warning('No precomputed data is found, not enough input arguments to compute data');
        return
    end
end

if exist(fname_out,'file') && force == false
    fprintf('Loading results template detection: %s\n', fname_out);
    % repeat to deal with load errors
    count = 0;
    err_count = 0;
    while count == err_count
        try
            load(fname_out, 'MuseStruct', 'Tindx_unique', 'LFP_avg');
        catch ME
            err_count = err_count + 1;
        end
        count = count + 1;
    end
    return
end

fprintf('Detecting spikes\n');

% get file format
[isNeuralynx, isMicromed, isBrainvision] = get_data_format(cfg);

% process templates
for itemp = 1 : size(template, 2)
    cfgtemp         = [];
    cfgtemp.latency = cfg.template.latency;
    template{itemp} = ft_selectdata(cfgtemp, template{itemp});

    if isfield(cfg.template,'resamplefs')
        cfgtemp             = [];
        cfgtemp.resamplefs  = cfg.template.resamplefs;
        template{itemp}     = ft_resampledata(cfgtemp, template{itemp});
    end
end

% loop over parts
for ipart = 1 :  size(cfg.directorylist,2)

    offset = 0; % in samples
    dirindx = []; %
    datlength = [];

    % loop over directories
    for idir = 1 : size(cfg.directorylist{ipart}, 2)

        if isNeuralynx
            nfile = size(cfg.cluster.channel, 2); % one file per channel
        elseif isMicromed
            nfile = 1; % only one file with all electrodes
            fname = fullfile(cfg.rawdir, [cfg.directorylist{ipart}{idir} '.TRC']);
        elseif isBrainvision
            nfile = 1; % only one file with all electrodes
            fname = fullfile(cfg.rawdir, [cfg.directorylist{ipart}{idir} '.eeg']);
        end

        % concatinate channels
        for ifile = 1 : nfile
            cfgtemp             = [];
            
            if isNeuralynx
                temp                = dir(fullfile(cfg.rawdir, cfg.directorylist{ipart}{idir}, ['*', cfg.cluster.channel{ifile},'.ncs']));
                cfgtemp.dataset     = fullfile(cfg.rawdir, cfg.directorylist{ipart}{idir}, temp.name);
                filedat{ifile}      = ft_preprocessing(cfgtemp);
                
                % labels have to be the same to append over directories
                temp                    = strsplit(filedat{ifile}.label{1},'_');
                filedat{ifile}.label{1} = strcat(temp{end-1},'_',temp{end});
            else
                cfgtemp.dataset     = fname;
                cfgtemp.channel     = cfg.cluster.channel;
                filedat{ifile}      = ft_preprocessing(cfgtemp);
            end
        end
        
        cfgtemp                     = [];
        cfgtemp.keepsampleinfo      = 'no';
        dirdat{idir}                = ft_appenddata(cfgtemp,filedat{:});
        clear filedat

        % bipolar rereferencing if requested
        % TODO: create template from scratch using indices from clusterLFP
        % then do same re-referencing as data
        if strcmp(cfg.template.reref, 'yes')
            if strcmp(cfg.template.refmethod, 'bipolar')
                labels_nonum    = regexprep(dirdat{idir}.label, '[0-9_]', '');
                [~,~,indx]      = unique(labels_nonum);
                clear group
                for i = 1 : max(indx)
                    cfgtemp             = [];
                    cfgtemp.reref       = 'yes';
                    cfgtemp.refmethod   = 'bipolar';
                    cfgtemp.channel     = dirdat{idir}.label(indx==i);
                    group{i}            = ft_preprocessing(cfgtemp,dirdat{idir});
                end
                dirdat{idir} = ft_appenddata([],group{:});
                clear group
            end
        end

        if isfield(cfg.template, 'resamplefs')
            cfgtemp            = [];
            cfgtemp.resamplefs = cfg.template.resamplefs;
            dirdat{idir}       = ft_resampledata(cfgtemp, dirdat{idir});
        end

        dirindx                = [dirindx ones(1,length(dirdat{idir}.time{1})) * idir];
        datlength(idir)        = length(dirdat{idir}.time{1});

        if idir > 1
            offset             = offset + length(dirdat{idir-1}.time{1});
            cfgtemp            = [];
            cfgtemp.offset     = offset;
            dirdat{idir}       = ft_redefinetrial(cfgtemp, dirdat{idir});
        end
    end

    cumsumdatlength = cumsum(datlength);
    cumsumdatlength = [0 cumsumdatlength(1:end-1)];

    cfgtemp                     = [];
    cfgtemp.keepsampleinfo      = 'no';
    temp                        = ft_appenddata(cfgtemp, dirdat{:});
    clear dirdat

    dat                         = rmfield(temp,'cfg');
    dat.trial                   = {cat(2, temp.trial{:})};
    dat.time                    = {cat(2, temp.time{:})};
    clear temp

    % sometimes data contains nan's probably at edges
    dat.trial{1}(isnan(dat.trial{1})) = 0;

    % find templates
    for itemp = 1 : size(template, 2)
        C                                           = normxcorr2(template{itemp}.avg', dat.trial{1}');
        C(1:size(template{itemp}.avg,2),:)          = nan;
        C(end-size(template{itemp}.avg,2):end,:)    = nan;
        noYshift                                    = size(template{itemp}.avg,1);
        C_norm{itemp}                               = normalize(C(:, noYshift)); clear C
        threshold                                   = nanstd(C_norm{itemp}) * cfg.template.threshold;
        [~, Tindx{ipart}{itemp}, ~, ~]              = findpeaks(C_norm{itemp}, 'MinPeakHeight', threshold, 'MinPeakDistance', dat.fsample*(-template{itemp}.time(1)+template{itemp}.time(end))/2);
    end

   % remove overlapping templates, selecting highest r
    a       = cell2mat(C_norm);
    a       = reshape(a, 1, numel(a));
    a       = normalize(a);
    a       = reshape(a, size(C_norm{1}, 1), size(C_norm, 2));
    C_norm2 = [];
    for i = 1 : size(C_norm, 2)
        C_norm2{i} = a(:, i);
    end
    clear a
    
    Tindx_unique{ipart} = Tindx{ipart};
    for itemp = 1 : size(template, 2)
        others          = 1 : size(template, 2);
        others(itemp)   = [];
        for iother = others
            [indxA, indxB]  = CommonElemTol(Tindx_unique{ipart}{itemp}, Tindx_unique{ipart}{iother}, 25);
            comp            = C_norm2{itemp}(indxA) >= C_norm2{iother}(indxB);
            Tindx_unique{ipart}{itemp}(indxA(~comp)) = [];
            Tindx_unique{ipart}{iother}(indxB(comp)) = [];
        end
    end

    % plot correlation, threshold and average LFG
    for itemp = 1 : size(template, 2)

        % plot threshold per template
        fig = figure('visible', cfg.visible); hold;
        plot(C_norm{itemp});
        axis tight
        ax = axis;
        plot([ax(1),ax(2)],[threshold, threshold],':k');
        if ~isempty(Tindx_unique{ipart}{itemp})
            scatter3(Tindx_unique{ipart}{itemp}, C_norm{itemp}(Tindx_unique{ipart}{itemp}), ones(size(Tindx_unique{ipart}{itemp}))*10, 'r.');
            n = size(Tindx_unique{ipart}{itemp}, 1);
        else
            n = 0;
        end
        axis tight
        box off
        title(sprintf('n = %d', n));

        % show separation in files and time
        ax = axis;
        axisnames = [];
        for i = 1 : length(cumsumdatlength)
            plot3([cumsumdatlength(i), cumsumdatlength(i)], [ax(3), ax(4)], [5, 5], 'color',[0, 0, 0]);
            axisnames{i} = datestr(MuseStruct{ipart}{i}.starttime);
        end
        set(gca, 'TickDir', 'out');
        xticks(cumsumdatlength);
        xticklabels(axisnames);
        xtickangle(90);

        % print to file
        set(fig,'PaperOrientation','landscape');
        set(fig,'PaperUnits','normalized');
        set(fig,'PaperPosition', [0 0 1 1]);
        fname_fig = fullfile(cfg.imagesavedir, 'templates', [cfg.prefix, 'p', num2str(ipart), '_template', num2str(itemp),'_threshold.png']);
        isdir_or_mkdir(fileparts(fname_fig));
        exportgraphics(fig, fname_fig);
        
        % skip further plotting if no templates were detected
        if isempty(Tindx_unique{ipart}{itemp})
            LFP_avg{ipart}{itemp} = [];
            continue
        end

        % create LFP averages
        startsample                     = Tindx_unique{ipart}{itemp} - size(template{itemp}.avg,2);
        endsample                       = Tindx_unique{ipart}{itemp};
        offset                          = zeros(size(Tindx_unique{ipart}{itemp})) + template{itemp}.time(1) * dat.fsample;
        cfgtemp                         = [];
        cfgtemp.trl                     = [startsample, endsample, offset];
        cfgtemp.trl                     = round(cfgtemp.trl);
        LFP_sel                         = ft_redefinetrial(cfgtemp, dat);
        LFP_avg{ipart}{itemp}           = ft_timelockanalysis([], LFP_sel);

        % plot LFPs vs. template
        fig = figure('visible', cfg.visible);

        subplot(1, 3, 1); hold;
        maxabs = -inf;
        for itrial = 1 : size(LFP_sel.trial,2)
            if max(max(abs(LFP_sel.trial{itrial}))) > maxabs
                maxabs = max(max(abs(LFP_sel.trial{itrial})));
            end
        end
        ytick = [];
        for itrial = 1 : size(LFP_sel.trial,2)
            ytick = [ytick, itrial * maxabs];
            for ichan = 1 : size(LFP_sel.label,1)
                plot(LFP_sel.time{itrial}, LFP_sel.trial{itrial}(ichan,:) + maxabs * ichan, 'color', [0.6, 0.6, 0.6]);
            end
        end
        for ichan = 1 : size(LFP_sel.label,1)
            plot(LFP_avg{ipart}{itemp}.time, LFP_avg{ipart}{itemp}.avg(ichan,:) + maxabs * ichan, 'color', [0, 0, 0]);
        end
        title(sprintf('n = %d', size(LFP_sel.trial,2)));
        yticks(ytick);
        set(gca,'TickLabelInterpreter', 'none')
        set(gca,'fontsize', 6)
        yticklabels(LFP_sel.label);
        axis tight
        box off

        % plot overlapping trial LFPs
        subplot(1, 3, 2); hold;
        maxabs = -inf;
        for itrial = 1 : size(LFP_avg{ipart}{itemp}.label, 1)
            if max(max(abs(LFP_avg{ipart}{itemp}.avg))) > maxabs
                maxabs = max(max(abs(LFP_avg{ipart}{itemp}.avg)));
            end
        end
        for ichan = 1 : size(LFP_avg{ipart}{itemp}.label,1)
            plot(template{itemp}.time, template{itemp}.avg(ichan,:) + maxabs * ichan, 'r');
            plot(LFP_avg{ipart}{itemp}.time, LFP_avg{ipart}{itemp}.avg(ichan,:) + maxabs * ichan, 'k');
        end
        set(gca,'fontsize', 6)
        set(gca,'Yticklabels',[])
        yticklabels([]);
        axis tight
        box off

        % plot average LPF template
        subplot(2, 3, 3);
        plot(template{itemp}.time, template{itemp}.avg', 'r');
        title('Template');
        set(gca,'fontsize', 6)
        axis tight
        box off

        % plot average LPF matches
        subplot(2, 3, 6);
        plot(LFP_avg{ipart}{itemp}.time, LFP_avg{ipart}{itemp}.avg', 'k');
        title('Average');
        set(gca,'fontsize', 6)
        axis tight
        box off

        % print to file
        set(fig,'PaperOrientation','landscape');
        set(fig,'PaperUnits','normalized');
        set(fig,'PaperPosition', [0 0 1 1]);
        fname_fig = fullfile(cfg.imagesavedir, 'templates', [cfg.prefix, 'p', num2str(ipart), '_template', num2str(itemp),'_LFP.png']);
        isdir_or_mkdir(fileparts(fname_fig));
        exportgraphics(fig, fname_fig);
        
    end

    % plot templates and thresholds for all templates
    fig = figure('visible', cfg.visible);

    for itemp = 1 : size(C_norm2, 2)

        subplot(size(C_norm2, 2) + 1, 1, itemp); hold;

        plot(C_norm2{itemp}, 'color', [0.5, 0.5, 0.5]);
        axis tight
        ax = axis;
        plot([ax(1),ax(2)],[threshold, threshold],':k');
        if ~isempty(Tindx_unique{ipart}{itemp})
            scatter3(Tindx{ipart}{itemp},           C_norm2{itemp}(Tindx{ipart}{itemp}),         ones(size(Tindx{ipart}{itemp}))*10, 'k.');
            scatter3(Tindx_unique{ipart}{itemp},    C_norm2{itemp}(Tindx_unique{ipart}{itemp}),  ones(size(Tindx_unique{ipart}{itemp}))*10, 'r.');
            n = size(Tindx_unique{ipart}{itemp}, 1);
        else
            n = 0;
        end
        title(sprintf('template%d, n = %d', itemp, n));
        set(gca, 'Layer', 'top');

        xticks('');
        set(gca,'fontsize', 6)
        for i = 1 : length(cumsumdatlength)
            plot3([cumsumdatlength(i), cumsumdatlength(i)], [ax(3), ax(4)], [5, 5], 'color', [0, 0, 0]);
            axisnames{i} = datestr(MuseStruct{ipart}{i}.starttime);
        end
        axis tight
        ax = axis;
        box off
    end

    % add filenames to extra subplot
    subplot(size(C_norm2, 2) + 1, 1, size(C_norm2, 2) + 1); hold;
    set(gca,'fontsize', 6)

    for i = 1 : length(cumsumdatlength)
        plot3([cumsumdatlength(i), cumsumdatlength(i)], [ax(3), ax(4)], [5, 5], 'color', [0, 0, 0]);
    end
    axis tight
    xlim([ax(1), ax(2)]);
    set(gca, 'TickDir', 'out');
    xticks(cumsumdatlength);
    set(gca, 'xticklabels', axisnames);
    xtickangle(90);

    % print to file
    set(fig,'PaperOrientation','landscape');
    set(fig,'PaperUnits','normalized');
    set(fig,'PaperPosition', [0 0 1 1]);    
    fname_fig = fullfile(cfg.imagesavedir, 'templates', [cfg.prefix, 'p', num2str(ipart), '_all_templates_threshold.png']);
    isdir_or_mkdir(fileparts(fname_fig));
    exportgraphics(fig, fname_fig);
    
    % add to MuseStruct and add to markerfile if requested
    for idir = unique(dirindx)

        % remove previous template markers
        names = fieldnames(MuseStruct{ipart}{idir}.markers);
        for itemp = 1 : 100
            for iname = 1 : size(names,1)
                if strcmp(names(iname), sprintf('template%d', itemp))
                    MuseStruct{ipart}{idir}.markers = rmfield(MuseStruct{ipart}{idir}.markers, sprintf('template%d', itemp));
                end
            end
        end

        % add new template markers
        for itemp = 1 : size(template, 2)

            % remove offset of directory
            indx = Tindx_unique{ipart}{itemp}(dirindx(Tindx_unique{ipart}{itemp}) == idir) - cumsumdatlength(idir) ;

            % remove offset of xcorr
            indx = indx - size(template{itemp}.time, 2);

            % shift to t=0 in template
            indx = indx - template{itemp}.time(1) * dat.fsample;

            % add to musestruct
            name = sprintf('template%d', itemp);
            MuseStruct{ipart}{idir}.markers.(name).comment       = 'Added by detectTemplate (Stephen)';
            MuseStruct{ipart}{idir}.markers.(name).color         = 'green';
            MuseStruct{ipart}{idir}.markers.(name).offset        = 'Added by detectTemplate (Stephen)';
            MuseStruct{ipart}{idir}.markers.(name).classgroupid  = '+3';
            MuseStruct{ipart}{idir}.markers.(name).editable      = 'Yes';
            MuseStruct{ipart}{idir}.markers.(name).classid       = '+666'; % will be replaced by writeMuseMarkers.m
            MuseStruct{ipart}{idir}.markers.(name).synctime      = (indx / dat.fsample)';
            MuseStruct{ipart}{idir}.markers.(name).clock         = (seconds(indx / dat.fsample) + MuseStruct{ipart}{idir}.starttime)';
        end

        if cfg.template.writemuse == true

            fname_mrk = fullfile(cfg.rawdir, cfg.directorylist{ipart}{idir},'Events.mrk');

            % backup markerfile
            if ~exist(cfg.muse.backupdir,'DIR')
                error('Backup directory does not exist');
            end
            [~, d] = fileparts(cfg.directorylist{ipart}{idir});
            if ~exist(fullfile(cfg.muse.backupdir, d), 'DIR')
                fprintf('Creating directory: %s\n', fullfile(cfg.muse.backupdir, d));
                eval(sprintf('!mkdir %s', fullfile(cfg.muse.backupdir, d)));
            end
            fname_backup = sprintf('Events_%s.mrk', datestr(now, 'mm-dd-yyyy_HH-MM-SS'));
            eval(sprintf('!cp %s %s', fname_mrk, fullfile(cfg.muse.backupdir, d, fname_backup)));
            fprintf('Succesfully backed up markerfile to %s\n',fullfile(cfg.muse.backupdir, d, fname_backup));

            % write to muse marker file
            writeMuseMarkerfile(MuseStruct{ipart}{idir}, fname_mrk);
        end
    end
    
    % regularly update inside loop, just in case
    save(fname_out, 'MuseStruct', 'Tindx_unique', 'LFP_avg', '-v7.3');
    close all
    clear LFP_sel C_norm C_norm2
end
