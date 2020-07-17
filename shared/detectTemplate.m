function [MuseStruct, C_norm, Tindx, LFP_avg] = detectTemplate(cfg, MuseStruct, template, force)

fname_out = fullfile(cfg.datasavedir,[cfg.prefix,'MuseStruct_detectedTemplates.mat']);
cfg.template.reref       = ft_getopt(cfg.template,'reref','no');
cfg.template.refmethod   = ft_getopt(cfg.template,'refmethod','none');
cfg.template.latency     = ft_getopt(cfg.template,'latency','all');
cfg.template.writemuse   = ft_getopt(cfg.template,'writemuse',true);
cfg.template.threshold   = ft_getopt(cfg.template,'threshold',3);
cfg.template.name        = ft_getopt(cfg.template,'name','TemplateDetect');
cfg.template.visible     = ft_getopt(cfg.template,'visible','on');

if exist(fname_out,'file') && force == false
    fprintf('************************************\n');
    fprintf('** Loading results spikedetection **\n');
    fprintf('************************************\n\n');
    load(fname_out, 'MuseStruct', 'C_norm', 'Tindx', 'LFP_avg');
    return
end


fprintf('**********************\n');
fprintf('** Detecting spikes **\n');
fprintf('**********************\n\n');

% get file format
[isNeuralynx, isMicromed, isBrainvision] = get_data_format(cfg);

% process template
cfgtemp         = [];
cfgtemp.latency = cfg.template.latency;
template        = ft_selectdata(cfgtemp,template);

if isfield(cfg.template,'resamplefs')
    cfgtemp             = [];
    cfgtemp.resamplefs  = cfg.template.resamplefs;
    template            = ft_resampledata(cfgtemp, template);
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
            if isNeuralynx
                temp                = dir(fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir},['*',cfg.cluster.channel{ifile},'.ncs']));
                fname{1}            = fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir},temp.name);
                filedat{ifile}      = ft_read_neuralynx_interp(fname);
                
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
        
        cfgtemp         = [];
        cfgtemp.latency = cfg.template.latency;
        template        = ft_selectdata(cfgtemp,template);
        
        if isfield(cfg.template,'resamplefs')
            cfgtemp            = [];
            cfgtemp.resamplefs = cfg.template.resamplefs;
            template           = ft_resampledata(cfgtemp, template);
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

    cfgtemp                     = [];
    cfgtemp.keepsampleinfo      = 'no';
    temp                        = ft_appenddata(cfgtemp, dirdat{:});
    clear dirdat

    dat                         = rmfield(temp,'cfg');
    dat.trial                   = {cat(2, temp.trial{:})};
    dat.time                    = {cat(2, temp.time{:})};
    clear temp
     
    % loop over templates
    
    % fitting
    fsample_dat_orig            = dat.fsample;
    fsample_template_orig       = round(1/mode(diff(template.time)));
%     C{ipart}                    = normxcorr2e(template.avg',dat.trial{1}','same');
    
    % sometimes data contains nan's probably at edges
    dat.trial{1}(isnan(dat.trial{1})) = 0; 
    C{ipart}                    = normxcorr2(template.avg', dat.trial{1}');

    % remove confounded beginning and end
    C{ipart}(1:size(template.avg,2),:)        = nan;
    C{ipart}(end-size(template.avg,2):end,:)  = nan;

    % plot detection points
    fig = figure('visible',cfg.template.visible); 
    for i = 1 : 7
        subplot(7,1,i);
        plot((C{ipart}(:,i)));
    end
    set(fig,'PaperOrientation','landscape');
    set(fig,'PaperUnits','normalized');
    set(fig,'PaperPosition', [0 0 1 1]);
    print(fig, '-dpng', fullfile(cfg.imagesavedir, [cfg.prefix, 'p', num2str(ipart), '_', cfg.template.name, '_detectTemplate_Xcorr.png']));
    print(fig, '-dpdf', fullfile(cfg.imagesavedir, [cfg.prefix, 'p', num2str(ipart), '_', cfg.template.name, '_detectTemplate_Xcorr.pdf']));    
    close all
    
    noYshift                = size(template.avg,1);
    C_norm{ipart}           = normalize(C{ipart}(:, noYshift));
    threshold               = nanstd(C_norm{ipart}) * cfg.template.threshold;
    [~,Tindx{ipart},~,~]    = findpeaks(C_norm{ipart}, 'MinPeakHeight', threshold, 'MinPeakDistance', dat.fsample*(-template.time(1)+template.time(end))/2);
        
    % greate LFP averages
    startsample             = Tindx{ipart} - size(template.avg,2);
    endsample               = Tindx{ipart};
    offset                  = zeros(size(Tindx{ipart})) + template.time(1) * dat.fsample;
    cfgtemp                 = [];
    cfgtemp.trl             = [startsample, endsample, offset];
    cfgtemp.trl             = round(cfgtemp.trl);
    LFP_sel                 = ft_redefinetrial(cfgtemp, dat);
    LFP_avg{ipart}          = ft_timelockanalysis([],LFP_sel);
    
    % plot LFPs vs. template
    fig = figure('visible',cfg.template.visible); 
    
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
        plot(LFP_avg{ipart}.time, LFP_avg{ipart}.avg(ichan,:) + maxabs * ichan, 'color', [0, 0, 0]);
    end
    title(sprintf('n = %d', size(LFP_sel.trial,2)));
    yticks(ytick);
    set(gca,'TickLabelInterpreter', 'none')
    set(gca,'fontsize', 6)
    yticklabels(LFP_sel.label);
    axis tight
    
    subplot(1, 3, 2); hold;
    maxabs = -inf;
    for itrial = 1 : size(LFP_avg{ipart}.label, 1)
        if max(max(abs(LFP_avg{ipart}.avg))) > maxabs
            maxabs = max(max(abs(LFP_avg{ipart}.avg)));
        end
    end
    for ichan = 1 : size(LFP_avg{ipart}.label,1)
        plot(template.time, template.avg(ichan,:) + maxabs * ichan,'r');
        plot(LFP_avg{ipart}.time, LFP_avg{ipart}.avg(ichan,:) + maxabs * ichan,'k');
    end
    set(gca,'TickLabelInterpreter', 'none')
    set(gca,'fontsize', 6)
    yticklabels(LFP_sel.label);
    axis tight
    box off
    
    subplot(2, 3, 3); 
    plot(template.time, template.avg','r');
    title('Template');
    subplot(2, 3, 6); 
    plot(LFP_avg{ipart}.time, LFP_avg{ipart}.avg','k');
    title('Average');

    set(fig,'PaperOrientation','landscape');
    set(fig,'PaperUnits','normalized');
    set(fig,'PaperPosition', [0 0 1 1]);
    print(fig, '-dpng', fullfile(cfg.imagesavedir, [cfg.prefix, 'p', num2str(ipart), '_', cfg.template.name, '_detectTemplate_LFP.png']));
    print(fig, '-dpdf', fullfile(cfg.imagesavedir, [cfg.prefix, 'p', num2str(ipart), '_', cfg.template.name, '_detectTemplate_LFP.pdf']));    
    close all
    
    fig = figure('visible',cfg.template.visible); hold;
    plot(C_norm{ipart});
    axis tight
    ax = axis;
    plot([ax(1),ax(2)],[threshold, threshold],':k');
    scatter3(Tindx{ipart}, C_norm{ipart}(Tindx{ipart}), ones(size(Tindx{ipart}))*10,'r.');
    axis tight
    box off
    
    cumsumdatlength = [0 cumsumdatlength(1:end-1)];
    
    % show separation in files and time
    for i = 1 : length(cumsumdatlength)
        plot3([cumsumdatlength(i), cumsumdatlength(i)], [ax(3), ax(4)], [-1, -1],'color',[0.7, 0.7, 0.7]);
        text(cumsumdatlength(i), ax(4), datestr(MuseStruct{ipart}{i}.starttime), 'rotation', 90, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'color', [0.7, 0.7, 0.7], 'fontsize', 8);
    end

    set(fig,'PaperOrientation','landscape');
    set(fig,'PaperUnits','normalized');
    set(fig,'PaperPosition', [0 0 1 1]);
    print(fig, '-dpng', fullfile(cfg.imagesavedir, [cfg.prefix, 'p', num2str(ipart), '_', cfg.template.name, '_detectTemplate_threshold.png']));
    print(fig, '-dpdf', fullfile(cfg.imagesavedir, [cfg.prefix, 'p', num2str(ipart), '_', cfg.template.name, '_detectTemplate_threshold.pdf']));   
    close all

    % indexes reflect end of template window
%     Tindx{ipart} = Tindx{ipart} - size(template.avg, 2) - template.time(1) * cfg.template.resamplefs;
    
    % add to MuseStruct and write    
    for idir = unique(dirindx)
        
        %  read Muse events file
        name_mrk    = fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir},'Events.mrk');
        
        % backup markerfile
        if ~exist(cfg.muse.backupdir,'DIR')
            error('Backup directory does not exist');
        end
        [~, d] = fileparts(cfg.directorylist{ipart}{idir});
        if ~exist(fullfile(cfg.muse.backupdir, d), 'DIR')
            fprintf('Creating directory: %s\n', fullfile(cfg.muse.backupdir,d));
            eval(sprintf('!mkdir %s', fullfile(cfg.muse.backupdir,d)));
        end
        fname_backup = sprintf('Events_%s.mrk', datestr(now, 'mm-dd-yyyy_HH-MM-SS'));
        eval(sprintf('!cp %s %s', name_mrk, fullfile(cfg.muse.backupdir, d, fname_backup)));
        fprintf('Succesfully backed up markerfile to %s\n',fullfile(cfg.muse.backupdir, d, fname_backup));
        
        % remove offset of directory     
        indx = Tindx{ipart}(dirindx(Tindx{ipart}) == idir) - cumsumdatlength(idir) ;
        
        % remove offset of xcorr 
        indx = indx - size(template.time,2);
        
        % shift to t=0 in template
        indx = indx - template.time(1) * dat.fsample;
        
        % add to musestruct
        MuseStruct{ipart}{idir}.markers.(cfg.template.name).comment       = 'Added by detectTemplate (Stephen)';
        MuseStruct{ipart}{idir}.markers.(cfg.template.name).color         = 'green';
        MuseStruct{ipart}{idir}.markers.(cfg.template.name).offset        = 'Added by detectTemplate (Stephen)';
        MuseStruct{ipart}{idir}.markers.(cfg.template.name).classgroupid  = '+3';
        MuseStruct{ipart}{idir}.markers.(cfg.template.name).editable      = 'Yes';
        MuseStruct{ipart}{idir}.markers.(cfg.template.name).classid       = '+666'; % will be replaced by writeMuseMarkers.m
        MuseStruct{ipart}{idir}.markers.(cfg.template.name).synctime      = (indx / dat.fsample)';
        MuseStruct{ipart}{idir}.markers.(cfg.template.name).clock         = (seconds(indx / dat.fsample) + MuseStruct{ipart}{idir}.starttime)';

        if cfg.template.writemuse == true
            writeMuseMarkerfile(MuseStruct{ipart}{idir},name_mrk);
        end
    end
    clear LFP_sel
end

save(fname_out,'MuseStruct', 'C_norm', 'Tindx', 'LFP_avg', '-v7.3');
