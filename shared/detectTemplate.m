function C = detectTemplate(cfg,MuseStruct,template,force)

fname_out = fullfile(cfg.datasavedir,[cfg.prefix,'MuseStruct_detectedTemplates.mat']);

if exist(fname_out,'file') && force == false
    fprintf('************************************\n');
    fprintf('** Loading results spikedetection **\n');
    fprintf('************************************\n\n');
    load(fullfile(cfg.datasavedir,[cfg.prefix,'MuseStruct_detectedTemplates.mat']));
    continue
end

    if force == true
        fprintf('**************************************\n');
        fprintf('** Forced redoing of templatedetection **\n');
        fprintf('**************************************\n\n');
    else
        fprintf('**********************\n');
        fprintf('** Detecting spikes **\n');
        fprintf('**********************\n\n');
    end

    % get file format
    [isNeuralynx, isMicromed, isBrainvision] = get_data_format(cfg);

    % loop over parts within subject
    for ipart = 1 : size(cfg.directorylist,2)

        fig = figure;

        % loop over directories
        for idir = 1 : size(cfg.directorylist{ipart}, 2)
            
            if isNeuralynx
                nfile = size(cfg.cluster.channel,2); % one file per channel
            elseif isMicromed
                nfile = 1; % only one file with all electrodes
                fname = fullfile(cfg.rawdir,[cfg.directorylist{ipart}{idir} '.TRC']);
            elseif isBrainvision
                nfile = 1; % only one file with all electrodes
                fname = fullfile(cfg.rawdir,[cfg.directorylist{ipart}{idir} '.eeg']);
            end

            % loop over files
            for ifile = 1 : nfile

                %load data
                if isNeuralynx
                    temp                = dir(fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir},['*',cfg.cluster.channel{ifile},'.ncs']));
                    fname{1}            = fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir},temp.name);
                    filedat{ifile}      = ft_read_neuralynx_interp(fname);
                else
                    cfgtemp.dataset     = fname;
                    cfgtemp.channel     = cfg.cluster.channel;
                    filedat{ifile}      = ft_preprocessing(cfgtemp);
                end

            end

            % concatinate channels
            cfgtemp                     = [];
            cfgtemp.keepsampleinfo      = 'no';
            dat                         = ft_appenddata(cfgtemp,filedat{:});
            clear filedat*

            % bipolar rereferencing if requested
%             cfg.cluster.reref       = ft_getopt(cfg.cluster,'reref','no');
%             cfg.cluster.refmethod   = ft_getopt(cfg.cluster,'refmethod','none');
%             if strcmp(cfg.cluster.refmethod, 'bipolar')
%                 labels_nonum    = regexprep(dat.label, '[0-9_]', '');
%                 [~,~,indx]      = unique(labels_nonum);
%                 clear group
%                 for i = 1 : max(indx)
%                     cfgtemp             = [];
%                     cfgtemp.reref       = 'yes';
%                     cfgtemp.refmethod   = 'bipolar';
%                     cfgtemp.channel     = dat.label(indx==i);
%                     group{i}            = ft_preprocessing(cfgtemp,dat);
%                 end
%                 dat = ft_appenddata([],group{:});
%                 clear group
%             end
%
            template_ds = template;

            cfgtemp = [];
            cfgtemp.latency = [-0.2 0.5];
            template_ds = ft_selectdata(cfgtemp,template);

%             cfgtemp             = [];
%             cfgtemp.resamplefs  = 250;
%             template_ds         = ft_resampledata(cfgtemp, template);
%             dat                 = ft_resampledata(cfgtemp, dat);

            C{ipart}{idir} = normxcorr2e(template_ds.avg',dat.trial{1}','same');
%             C{ipart}{idir} = xcorr2(dat.trial{1}',template.trial{1}');


            noYshift            = floor(size(C{ipart}{idir},2) / 2) + 1;
            C{ipart}{idir}      = C{ipart}{idir}(:, noYshift);
            C{ipart}{idir}      = normalize(C{ipart}{idir});
            threshold           = std(C{ipart}{idir})*3;
            [Ypk,Xpk,Wpk,Ppk]   = findpeaks(C{ipart}{idir},'MinPeakHeight', threshold, 'MinPeakDistance', dat.fsample*(-template.time(1)+template.time(end))/2);

            cfgtemp = [];
            cfgtemp.trl = [Xpk+template.time(1)*dat.fsample, Xpk+template.time(end)*dat.fsample, zeros(size(Xpk))+template.time(1)*dat.fsample];
            LFP_sel = ft_redefinetrial(cfgtemp, dat);

            maxabs = -inf;
            for itrial = 1 : size(LFP_sel.trial,2)
                if max(max(abs(LFP_sel.trial{itrial}))) > maxabs
                    maxabs = max(max(abs(LFP_sel.trial{itrial})));
                end
            end

            subplot(1,size(cfg.directorylist{ipart}, 2),idir);
            hold;
            for itrial = 1 : size(LFP_sel.trial,2)
                plot(LFP_sel.time{itrial}, LFP_sel.trial{itrial}+maxabs * itrial);
            end

            LFP_avg = ft_timelockanalysis([],LFP_sel);
%
%             figure;
%             subplot(2,1,1);
%             plot(template.time, template.avg);
%             subplot(2,1,2);
%             plot(LFP_avg.time, LFP_avg.avg')
%

        end
        set(fig,'PaperOrientation','landscape');
        set(fig,'PaperUnits','normalized');
        set(fig,'PaperPosition', [0 0 1 1]);
        print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,'p',num2str(ipart),'_detected_Template.png']));
        close all
    end
