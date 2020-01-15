function detectSpikes(cfg,MuseStruct_micro,MuseStruct_macro,force,writeMuseMarker)

fname_out = fullfile(cfg.datasavedir,[cfg.prefix,'MuseStruct_detectedSpikes.mat']);

if exist(fname_out,'file') && force == false
    fprintf('************************************\n');
    fprintf('** Loading results spikedetection **\n');
    fprintf('************************************\n\n');
    load(fullfile(cfg.datasavedir,[cfg.prefix,'MuseStruct_detectedSpikes.mat']));
else
    
    if force == true
        fprintf('**************************************\n');
        fprintf('** Forced redoing of spikedetection **\n');
        fprintf('**************************************\n\n');
    else
        fprintf('**********************\n');
        fprintf('** Detecting spikes **\n');
        fprintf('**********************\n\n');
    end
    
    for ipart = 1:length(MuseStruct_macro)   
        
        for idir = 1:length(MuseStruct_macro{ipart})       
            
            % select MICRO files
            micro_filenrs = [];
            for ifile = 1 : size(MuseStruct_micro{ipart}{idir}.filenames,2)
                for ilabel = 1 : size(cfg.labels.micro,2)
                    if ~isempty(strfind(MuseStruct_micro{ipart}{idir}.filenames{ifile},cfg.labels.micro{ilabel}))
                        micro_filenrs       = [micro_filenrs, ifile];
                        microlabel{ifile}   = cfg.labels.micro{ilabel};
                    end
                end
            end
            
            % select MACRO files
            macro_filenrs = [];
            for ifile = 1 : size(MuseStruct_macro{ipart}{idir}.filenames,2)
                for ilabel = 1 : size(cfg.labels.macro,2)
                    if ~isempty(strfind(MuseStruct_macro{ipart}{idir}.filenames{ifile},cfg.labels.macro{ilabel}))
                        macro_filenrs       = [macro_filenrs, ifile];
                        macrolabel{ifile}   = cfg.labels.macro{ilabel};
                    end
                end
            end
            
            hasdata_micro = true(size(micro_filenrs));
            hasdata_macro = true(size(macro_filenrs));
            %
            %             % load trials for selected MICRO channels
            %             for ifile = micro_filenrs
            %
            %                 % to deal with missing data
            %                 fname{1}                = fullfile(MuseStruct_micro{ipart}{idir}.directory, MuseStruct_micro{ipart}{idir}.filenames{ifile});
            %                 try
            %                     dat{ifile}                     = ft_read_neuralynx_interp(fname);
            %             dat_micro{ifile}.label{1} = dat_micro{ifile}.label{1}(end-5:end);
            %
            %                 catch
            %                 end
            %             end
            
            % load trials for selected MACRO channels
            for ifile = macro_filenrs
                
                % to deal with missing data
                fname{1}                = fullfile(MuseStruct_macro{ipart}{idir}.directory, MuseStruct_macro{ipart}{idir}.filenames{ifile});
                try
                    dat_macro{ifile} = ft_read_neuralynx_interp(fname);
                    dat_macro{ifile}.label{1} = dat_macro{ifile}.label{1}(end-5:end);
                catch
                    fprintf('problems with file %s\n',fname{1});
                    hasdata_micro(ifile) = false;
                end
            end
            
            
            % concatinate channels, separately for MICRO/MACRO
            cfgtemp                             = [];
            %             cfgtemp.keepsampleinfo              = 'no';
            %             dirdat_micro{idir}                  = ft_appenddata(cfgtemp,filedat_micro{micro_filenrs(hasdata_micro)});
            dat_macro_append{idir}                  = ft_appenddata(cfgtemp,dat_macro{macro_filenrs(hasdata_macro)});
            clear dat_macro
            
            % rereference
%             cfgtemp = [];
%             cfgtemp.reref = 'yes';
%             cfgtemp.refmethod = 'bipolar';
%             dat_macro_append{idir} = ft_preprocessing(cfgtemp,dat_macro_append{idir} );

            % flip data positive negative
            dat_macro_append{idir}.trial{1} = -dat_macro_append{idir}.trial{1};
            
            
            LS = 5;         % Left half-wave slope; default: 7
            RS = 5;         % Right half-wave slope; default: 7
            TAMP = 400;     % Total amplitude; default: 600
            LD = 1;        % Left half-wave duration; default = 10
            RD = 1;        % Right half-wave duration; default = 10
            STDCoeff = 4;   % Chebyshev inequality coefficient (distance from centre point or mean); default 4
            SCALE = 70;     % Scaling parameter
            BlockSize = 1;  % Data processing block size in minutes
            TroughSearch = 40; % distance in ms to search for a trough on each side of a detected peak
            DetThresholds = [ LS; RS; TAMP; LD; RD;];
            FilterSpec =  [20; 50; 1; 35;];
            
            [SpikeIndex, ChanId, SpikeFV] = mDetectSpike(dat_macro_append{idir}.trial{1}',dat_macro_append{idir}.fsample, ...
                cfg.spikedetect.BlockSize,cfg.spikedetect.SCALE,cfg.spikedetect.STDCoeff,cfg.spikedetect.DetThresholds,cfg.spikedetect.FilterSpec,cfg.spikedetect.TroughSearch);
            
            % iterate to remove detections that are too close together
            for i = 1 : 10
                toremove = diff(SpikeIndex / dat_macro_append{idir}.fsample) < 0.100;
                SpikeIndex = SpikeIndex(~toremove);
                ChanId = ChanId(~toremove);
                SpikeFV = SpikeFV(~toremove,:);
            end

            if ~isempty(SpikeIndex)
                cfgtemp                         = [];
                cfgtemp.trl(:,1)                = SpikeIndex - 0.5 * dat_macro_append{idir}.fsample;
                cfgtemp.trl(:,2)                = SpikeIndex + 0.5 * dat_macro_append{idir}.fsample;
                cfgtemp.trl(:,3)                = - 0.5 * dat_macro_append{idir}.fsample;
                cfgtemp.trl(:,4)                = ChanId;
                dat_timelocked{idir}            = ft_redefinetrial(cfgtemp,dat_macro_append{idir});
                clear dat_macro_append
                
                avg_timelocked{idir}            = ft_timelockanalysis([],dat_timelocked{idir});
                
                
                % plot
                fig             = figure;
                fig.Renderer    = 'Painters'; % Else pdf is saved to bitmap
                
                subplot(2,1,1);
                hold;
                h               = 1200;
                for itrial = 1 : size(dat_timelocked{idir}.trial,2)
                    plot(dat_timelocked{idir}.time{itrial},dat_timelocked{idir}.trial{itrial}(dat_timelocked{idir}.trialinfo(itrial,1),:) + itrial*h);
                end
                axis tight
                plot([0,0],ylim,':k');
                subplot(2,1,2);
                plot(avg_timelocked{idir}.time,avg_timelocked{idir}.avg);
                xlabel('Time (s)');
                axis tight
                
                % print to file
                set(fig,'PaperOrientation','landscape');
                set(fig,'PaperUnits','normalized');
                set(fig,'PaperPosition', [0 0 1 1]);
                print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'p',num2str(ipart),'_d',num2str(idir),'-',MuseStruct_macro{ipart}{idir}.filenames{ifile}(1:end-6),'.pdf']));
%                 print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,'p',num2str(ipart),'_d',num2str(idir),'-',MuseStruct_macro{ipart}{idir}.filenames{ifile}(1:end-6),'.png']),'-r300');
                close all
                
                clear avg_timelocked dat_timelocked
                
                % write markers to muse
                
                % update marker
                MuseStruct_macro{ipart}{idir}.markers.SpikeDetect.offset = 0;
                if writeMuseMarker
                    
                    %  read Muse events file
                    name_mrk    = fullfile(MuseStruct_macro{ipart}{idir}.directory,'Events.mrk');
                    
                    % backup markerfile
                    if ~exist(cfg.muse.backupdir,'DIR')
                        error('Backup directory does not exist');
                    end
                    [~, d] = fileparts(MuseStruct_macro{ipart}{idir}.directory);
                    if ~exist(fullfile(cfg.muse.backupdir,d),'DIR')
                        fprintf('Creating directory: %s\n',fullfile(cfg.muse.backupdir,d));
                        eval(sprintf('!mkdir %s',fullfile(cfg.muse.backupdir,d)));
                    end
                    fname_backup = sprintf('Events_%s.mrk', datestr(now, 'mm-dd-yyyy_HH-MM-SS'));
                    eval(sprintf('!cp %s %s',name_mrk,fullfile(cfg.muse.backupdir,d,fname_backup)));
                    fprintf('Succesfully backed up markerfile to %s\n',fullfile(cfg.muse.backupdir,d,fname_backup));
                    
                    MuseStruct_macro{ipart}{idir}.markers.SpikeDetect.comment       = 'Added by detectSpikes (Stephen)';
                    MuseStruct_macro{ipart}{idir}.markers.SpikeDetect.color         = 'red';
                    MuseStruct_macro{ipart}{idir}.markers.SpikeDetect.offset        = 'Added by detectSpikes (Stephen)';
                    MuseStruct_macro{ipart}{idir}.markers.SpikeDetect.classgroupid  = '+3';
                    MuseStruct_macro{ipart}{idir}.markers.SpikeDetect.editable      = 'Yes';
                    MuseStruct_macro{ipart}{idir}.markers.SpikeDetect.classid       = '+666'; % will be replaced by writeMuseMarkers.m
                    MuseStruct_macro{ipart}{idir}.markers.SpikeDetect.synctime      = SpikeIndex / MuseStruct_macro{ipart}{idir}.Fs;
                    writeMuseMarkers(MuseStruct_macro{ipart}{idir},name_mrk);
                end
            end % ~isempty(SpikeIndx)
            
        end % idir
    end % ipart
end

% if savedat
%     % fname_out
% end
