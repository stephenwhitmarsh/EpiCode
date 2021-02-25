function [SpikeRaw, SpikeTrials] = readSpykingCircus_allmarkers3(cfg,MuseStruct)

addpath /network/lustre/iss01/charpier/stephen.whitmarsh/WhitmarshEpilepsy
addpath /network/lustre/iss01/charpier/stephen.whitmarsh/fieldtrip/
addpath /network/lustre/iss01/charpier/stephen.whitmarsh/WhitmarshEpilepsy/mlib6/
ft_defaults

fname = fullfile(cfg.datasavedir,'spikedata_concatinated.mat');
if exist(fname,'file') && cfg.force == false
    load(fname,'SpikeRaw','SpikeTrials');
else
    
    
    %% use spiking circus output
    
    % find spiking-circus path
    temp = dir(fullfile(cfg.datasavedir,'SpykingCircus',['all_concatinated_',cfg.channel(1:end-2),'_*.result',cfg.suffix,'.hdf5']));
    fname_spikes = fullfile(temp.folder,temp.name);
    
    temp = dir(fullfile(cfg.datasavedir,'SpykingCircus',['all_concatinated_',cfg.channel(1:end-2),'_*.templates',cfg.suffix,'.hdf5']));
    fname_templates = fullfile(temp.folder,temp.name);
    
    % load spiking data
    if exist(fname_spikes,'file')
        fprintf('Loading spike data from: %s\n',fname_spikes);
        datinfo     = h5info(fname_spikes);
        temp        = dir(fullfile(cfg.datasavedir,['all_concatinated_',cfg.channel(1:end-2),'_*.ncs']));
        hdr         = ft_read_header(fullfile(temp(1).folder,temp(1).name)); % take the first file to extract the header of the data
        timestamps  = ft_read_data(fullfile(temp(1).folder,temp(1).name),'timestamp','true');  % take the first concatinated file to extract the timestamps
        
        % read spiketimes of clusters
        for i = 1 : size(datinfo.Groups,1)
            names(i) = string(datinfo.Groups(i).Name);
            if strfind(names(i),'spiketimes')
                spiketimes_indx = i;
            end
        end
        
        for i = 1 : size(datinfo.Groups(spiketimes_indx).Datasets,1) % number of templates
            SpikeRaw.label{i} = datinfo.Groups(spiketimes_indx).Datasets(i).Name;
            temp = strsplit(SpikeRaw.label{i},'_');
            clusternr(i) = str2double(temp{2});
        end
        
        for i = 1:numel(clusternr)
            % read spike timings (in seconds)
            datasetname = char(strcat('/spiketimes/',SpikeRaw.label{i}));
            SpikeRaw.sample{clusternr(i)+1} = h5read(fname_spikes,datasetname); % count from 1 instead of 0
            
            % read amplitudes
            datasetname = char(strcat('/amplitudes/',SpikeRaw.label{i}));
            SpikeRaw.amplitude{clusternr(i)+1} = h5read(fname_spikes,datasetname); % count from 1 instead of 0
            
            % map samplenrs onto timestamps
            SpikeRaw.timestamp{i} = timestamps(SpikeRaw.sample{i});
            %         SpikeRaw.timestamp{i} = SpikeRaw.sample{i} * hdr.TimeStampPerSample + double(hdr.FirstTimeStamp);
            
        end
        
        % load templates
        templates_size  = double(h5read(fname_templates, '/temp_shape'));
        N_e             = templates_size(2);
        N_t             = templates_size(1);
        temp_x          = double(h5read(fname_templates, '/temp_x') + 1);
        temp_y          = double(h5read(fname_templates, '/temp_y') + 1);
        temp_z          = double(h5read(fname_templates, '/temp_data'));
        templates       = sparse(temp_x, temp_y, temp_z, templates_size(1)*templates_size(2), templates_size(3));
        templates_size  = [templates_size(1) templates_size(2) templates_size(3)/2];
        
        for itemp = 1:numel(clusternr)
            template = full(reshape(templates(:, itemp), templates_size(2), templates_size(1)))';
            [~,i] = max(mean(abs(template),2));
            SpikeRaw.template(itemp,:,:) = template;
            SpikeRaw.template_maxchan(itemp) = i;
        end
        
        
        for ilabel = 1 : size(cfg.label,2)
            
            % clock time of each event
            clocktimes = [];
            for ifile = 1 : size(MuseStruct,2)
                if isfield(MuseStruct{ifile}.markers,cfg.startend{ilabel})
                    if isfield(MuseStruct{ifile}.markers.(cfg.startend{ilabel}),'clock')
                        clocktimes = [clocktimes, MuseStruct{ifile}.markers.(cfg.startend{ilabel}).clock];
                    end
                end
            end
            
            % create spiketrials timelocked to events
            cfgtemp                         = [];
            cfgtemp.trl                     = cfg.trialinfo_single{ilabel} + cfg.trialinfo_all(ilabel,1) - 1; % + cfg.prestim(ilabel) * hdr.Fs;
            cfgtemp.trl(:,3)                = -ones(size(cfgtemp.trl,1),1) * cfg.prestim(ilabel) * hdr.Fs;
            cfgtemp.trlunit                 = 'samples';
            cfgtemp.hdr                     = hdr;
            SpikeTrials{ilabel}             = ft_spike_maketrials(cfgtemp,SpikeRaw);
            
            % spike density function
            cfgtemp                         = [];
            cfgtemp.fsample                 = cfg.resamplefs;   % sample at 1000 hz
            cfgtemp.timwin                  = cfg.timwin{ilabel}; % [-0.0025 0.0025];
            sdf                             = ft_spikedensity(cfgtemp,SpikeTrials{ilabel});
            
            % ISI
            cfgtemp                         = [];
            cfgtemp.bins                    = [0:0.0005:0.1];   % use bins of 0.5 milliseconds
            cfgtemp.param                   = 'coeffvar';       % compute the coefficient of variation (sd/mn of isis)
            isih                            = ft_spike_isi(cfgtemp,SpikeTrials{ilabel});
            
            %% plot rasterplot for each template
            for itemp = 1 : size(SpikeRaw.label,2)
                
                % start figure
                fig = figure;
                
                % select template and create timeaxis for template
                tempsel = squeeze(SpikeRaw.template(itemp,SpikeRaw.template_maxchan(itemp),:));
                temptime = ((1:size(SpikeRaw.template,3))/hdr.Fs*1000)';
                
                % peak detection according to MATLAB - have to plot if first,
                % else a bug occurs because some handles get deleted
                subplot(5,4,8);
                [~,~,W,~] = findpeaks(tempsel,temptime,'NPeaks',1,'SortStr','descend');
                findpeaks(tempsel,temptime,'NPeaks',1,'SortStr','descend','Annotate','extents');
                legend('off');
                title(sprintf('%.2fms',W));
                set(gca, 'XTickLabel', '');
                axis tight
                
                % Spike density + raster
                hs0                     = subplot(1,4,1);
                pos0                    = get(hs0, 'Position');
                cfgtemp                 = [];
                cfgtemp.spikechannel    = itemp;
                cfgtemp.topplotsize     = 0.2;
                cfgtemp.topplotfunc     = 'line'; % plot as a line
                cfgtemp.latency         = [-cfg.prestim(ilabel), cfg.poststim(ilabel)];
                % cfgtemp.spikelength     = 1;
                % cfgtemp.linewidth       = 1;
                cfgtemp.errorbars       = 'sem'; % plot with the standard deviation
                cfgtemp.interactive     = 'yes'; % toggle off interactive mode
                ft_spike_plot_raster(cfgtemp,SpikeTrials{ilabel}, sdf)
                
                title('Spike density and rasterplot');
                set(hs0,'YaxisLocation','left');
                ylabel('')
                
                % plot histogram of spike frequencies
                hs3 = subplot(1,4,2);
                pos3 = get(hs3, 'Position');
                pos3 = pos3 .* [0.9 1 1 0.8]; % and decrease the height to match rasterplot
                sel = (SpikeTrials{ilabel}.time{itemp} >= -cfg.prestim(ilabel) & SpikeTrials{ilabel}.time{itemp} <= cfg.poststim(ilabel));
                clear f;
                for itrial = 1 : size(SpikeTrials{ilabel}.trialtime,1)
                    f(itrial) = sum(SpikeTrials{ilabel}.trial{itemp}(sel) == itrial);
                end
                nr_per_bin = 10;
                indx = ceil((1 : size(SpikeTrials{ilabel}.trialtime,1))/nr_per_bin);
                f2 = splitapply(@(x1)(mean(x1)),f,indx);
                bar(f2',1,'stacked','k');
                
                %             bar(f,1,'facecolor',[0 0 0]);
                view(90,90)
                set(hs3, 'Position', pos3, 'XTickLabel', '','xtick',[],'YaxisLocation','left');
                axis tight
                grid on
                grid minor
                title('Spikecount');
                
                % plot clock time
                hs1 = subplot(1,4,3);
                pos1 = get(hs1, 'Position');
                plot(clocktimes,'-k'); grid on;
                axis tight
                view(90,90)
                pos1 = pos1 .* [0.85 1 1 0.8]; % and decrease the height to match rasterplot
                set(hs1, 'Position', pos1, 'XTickLabel', '');
                ytickangle(45);
                set(hs1,'Fontsize',8);
                xt = yticks;
                
                % plot histogram of event frequencies
                hs2 = subplot(6,4,3);
                pos2 = [pos1(1) (pos1(2)+pos1(4)) pos1(3) pos0(4)*0.2];
                [histcount,edges,bin] = histcounts(clocktimes,xt);
                bar(histcount,1,'facecolor',[0 0 0]);
                set(hs2, 'Position', pos2, 'XTickLabel', '','xtick',[],'YaxisLocation','right');
                axis tight
                title('Eventcount');
                
                % template
                subplot(5,4,4); hold;
                maxheight = max(max(abs(SpikeRaw.template(itemp,:,:))));
                for ichan = 1 : size(SpikeRaw.template,2)
                    if ichan == SpikeRaw.template_maxchan(itemp)
                        plot(temptime,squeeze(SpikeRaw.template(itemp,ichan,:))+(ichan)*maxheight,'r');
                    else
                        plot(temptime,squeeze(SpikeRaw.template(itemp,ichan,:))+(ichan)*maxheight,'k');
                    end
                end
                axis tight
                xlabel('ms');
                ax = axis;
                axis([ax(1), ax(2), 0, ax(4)]);
                title('Template');
                set(gca, 'YTickLabel', '','XTickLabel', '');
                
                % peak width accoridng to Gast et. al
                subplot(5,4,12); hold;
                plot(temptime,tempsel,'k');
                axis tight
                [Ypos,Xpos] = findpeaks(tempsel',temptime,'NPeaks',1,'SortStr','descend');
                [Yneg,Xneg] = findpeaks(-tempsel',temptime,'NPeaks',2,'SortStr','descend','MinPeakDistance',0.5);
                plot([Xpos,Xneg(1)],[Yneg(1)*0.3, Yneg(1)*0.3],'-o','Color',[1 0 0],'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]);
                plot([Xpos,Xneg(2)],-[Yneg(2)*0.3, Yneg(2)*0.3],'-o','Color',[0 0 1],'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]);
                title(sprintf('%.2fms, %.2fms',abs(Xpos-Xneg(1)),abs(Xpos-Xneg(2))));
                set(gca, 'XTickLabel', '');
                set(gca, 'YTickLabel', '');
                
                subplot(5,4,16);
                [fwhm1,fwhm2] = mwave(tempsel,1/hdr.Fs,'plot');
                ax = axis;
                d = 5;
                set(gca, 'XTickLabel', '');
                set(gca, 'XTick', '');
                set(gca, 'YTickLabel', '');
                title(sprintf('%.2fms, %.2fms',fwhm2*1000,fwhm1*1000))
                ylabel(''); xlabel('');
                
                % ISI + ISI return plot
                subplot(5,4,20); hold;
                %         cfgtemp                 = [];
                %         cfgtemp.spikechannel    = isih.label{itemp};
                %         cfgtemp.interpolate     = 5; % interpolate at 5 times the original density
                %         cfgtemp.window          = 'gausswin'; % use a gaussian window to smooth
                %         cfgtemp.winlen          = 0.004; % the window by which we smooth has size 4 by 4 ms
                %         cfgtemp.colormap        = parula(300); % colormap
                %         cfgtemp.scatter         = 'no'; % do not plot the individual isis per spike as scatters
                %         ft_spike_plot_isireturn(cfgtemp,isih)
                bar(isih.time,isih.avg(itemp,:),1);
                [y,indx] = max(isih.avg(itemp,:));
                title(sprintf('Max: %.4fms',isih.time(indx)));
                xlabel('ms');
                axis tight
                grid on
                
                % write to PDF
                fig.Renderer='Painters';
                set(fig,'PaperOrientation','landscape');
                set(fig,'PaperUnits','normalized');
                set(fig,'PaperPosition', [0 0 1 1]);
                print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.label{ilabel},'_raster_template_fromall_',num2str(itemp),'.pdf']),'-r600');
            end
            
            %% plot overview of all clusters
            fig = figure;
            
            % ISI
            cfgtemp                         = [];
            cfgtemp.bins                    = [0:0.001:0.05];   % use bins of 0.5 milliseconds
            cfgtemp.param                   = 'coeffvar';       % compute the coefficient of variation (sd/mn of isis)
            isih                            = ft_spike_isi(cfgtemp,SpikeTrials{ilabel});
            
            for itemp = 1 : size(SpikeRaw.label,2)
                
                % select template and create timeaxis for template
                tempsel = squeeze(SpikeRaw.template(itemp,SpikeRaw.template_maxchan(itemp),:));
                temptime = ((1:size(SpikeRaw.template,3))/hdr.Fs*1000)';
                
                % template
                subplot(size(SpikeRaw.label,2),4,(itemp-1)*4+1); hold;          
                maxheight = max(max(abs(SpikeRaw.template(itemp,:,:))));
                for ichan = 1 : size(SpikeRaw.template,2)
                    if ichan == SpikeRaw.template_maxchan(itemp)
                        plot(temptime,squeeze(SpikeRaw.template(itemp,ichan,:))+(ichan)*maxheight,'r');
                    else
                        plot(temptime,squeeze(SpikeRaw.template(itemp,ichan,:))+(ichan)*maxheight,'k');
                    end
                end
                axis tight
                xlabel('ms');
                ax = axis;
                axis([ax(1), ax(2), 0, ax(4)]);
                grid on
                set(gca, 'YTickLabel', '');

                % peak width accoridng to Gast et. al
                subplot(size(SpikeRaw.label,2),4,(itemp-1)*4+2); hold;
                plot(temptime,tempsel,'k');
                axis tight
                [Ypos,Xpos] = findpeaks(tempsel',temptime,'NPeaks',1,'SortStr','descend');
                [Yneg,Xneg] = findpeaks(-tempsel',temptime,'NPeaks',2,'SortStr','descend','MinPeakDistance',0.5);
                plot([Xpos,Xneg(1)],[Yneg(1)*0.3, Yneg(1)*0.3],'-o','Color',[1 0 0],'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]);
                plot([Xpos,Xneg(2)],-[Yneg(2)*0.3, Yneg(2)*0.3],'-o','Color',[0 0 1],'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]);
                title(sprintf('%.2fms, %.2fms',abs(Xpos-Xneg(1)),abs(Xpos-Xneg(2))));
%                 set(gca, 'XTickLabel', '');
                set(gca, 'YTickLabel', '');
                
                subplot(size(SpikeRaw.label,2),4,(itemp-1)*4+3); hold;
                [fwhm1,fwhm2] = mwave(tempsel,1/hdr.Fs,'plot');
                ax = axis;
                d = 5;
                set(gca, 'XTickLabel', '');
                set(gca, 'XTick', '');
                set(gca, 'YTickLabel', '');
                title(sprintf('%.2fms, %.2fms',fwhm2*1000,fwhm1*1000))
                ylabel(''); xlabel('');
                
                % ISI + ISI return plot
                subplot(size(SpikeRaw.label,2),4,(itemp-1)*4+4); hold;
                %         cfgtemp                 = [];
                %         cfgtemp.spikechannel    = isih.label{itemp};
                %         cfgtemp.interpolate     = 5; % interpolate at 5 times the original density
                %         cfgtemp.window          = 'gausswin'; % use a gaussian window to smooth
                %         cfgtemp.winlen          = 0.004; % the window by which we smooth has size 4 by 4 ms
                %         cfgtemp.colormap        = parula(300); % colormap
                %         cfgtemp.scatter         = 'no'; % do not plot the individual isis per spike as scatters
                %         ft_spike_plot_isireturn(cfgtemp,isih)
                bar(isih.time,isih.avg(itemp,:),1);
                [y,indx] = max(isih.avg(itemp,:));
                title(sprintf('Max: %.4fms',isih.time(indx)));
                xlabel('ms');
                axis tight
                grid on
            end

            % write to PDF
            fig.Renderer='Painters';
            set(fig,'PaperOrientation','landscape');
            set(fig,'PaperUnits','normalized');
            set(fig,'PaperPosition', [0 0 1 1]);
            print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.label{ilabel},'_all_templates_morphology.pdf']));
            
            
            %% plot raster for all clusters
            fig = figure;
            
            % Spike density + raster
            cmap = lines(size(SpikeTrials{ilabel}.label,2));
            subplot(1,3,1);
            cfgtemp                 = [];
            cfgtemp.spikechannel    = 'all';
            cfgtemp.topplotsize     = 0.2;
            cfgtemp.topplotfunc     = 'line'; % plot as a line
            cfgtemp.latency         = [-cfg.prestim(ilabel), cfg.poststim(ilabel)];
            %     cfgtemp.spikelength     = 1;
            %     cfgtemp.linewidth       = 1;
            cfgtemp.errorbars       = 'sem'; % plot with the standard deviation
            cfgtemp.interactive     = 'yes'; % toggle off interactive mode
            cfgtemp.cmapneurons     = cmap;
            ft_spike_plot_raster(cfgtemp,SpikeTrials{ilabel}, sdf)
            title('Spike density and rasterplot');
            %
            %         % plot clock time
            %         hs = subplot(1,3,2);
            %         plot(clocktimes,'-k'); grid on;
            %         pos = get(hs, 'Position');
            %         axis tight
            %         view(90,90)
            %         pos2 = pos .* [0.9 1 1 0.8]; % and decrease the height to match rasterplot
            %         set(hs, 'Position', pos2, 'XTickLabel', '','xtick',[]);
            %         ytickangle(45);
            %         set(hs,'Fontsize',8);
            
            
            % plot histogram of spike frequencies
            
            hs3 = subplot(1,3,2); hold;
            pos3 = get(hs3, 'Position');
            pos3 = pos3 .* [0.9 1 1 0.8]; % and decrease the height to match rasterplot
            clear f
            
            for itemp = 1 : size(SpikeTrials{ilabel}.label,2)
                sel = (SpikeTrials{ilabel}.time{itemp} > -cfg.prestim(ilabel) & SpikeTrials{ilabel}.time{itemp} < cfg.poststim(ilabel));
                for itrial = 1 : size(SpikeTrials{ilabel}.trialtime,1)
                    f(itemp,itrial) = sum(SpikeTrials{ilabel}.trial{itemp}(sel) == itrial);
                end
            end
            
            nr_per_bin = 10;
            indx = ceil((1 : size(SpikeTrials{ilabel}.trialtime,1))/nr_per_bin);
            f2 = splitapply(@(x1)(mean(x1,2)),f,indx);
            h = bar(f2',1,'stacked');
            set(h,{'FaceColor'},num2cell(cmap,2));
            
            %         f = squeeze(nansum(sdf.trial(:,itemp,:)>0,3)); % is for whole trial with variable triallengths
            view(90,90)
            set(hs3, 'Position', pos3, 'XTickLabel', '','xtick',[],'YaxisLocation','left');
            axis tight
            grid on
            grid minor
            title('Spikecount');
            
            % template
            subplot(1,3,3); hold;
            maxheight = max(max(max(abs(SpikeRaw.template(:,:,:)))));
            for itemp = 1 : size(SpikeTrials{ilabel}.label,2)
                for ichan = 1 : size(SpikeRaw.template,2)
                    plot([1:size(SpikeRaw.template,3)]*hdr.Fs,squeeze(SpikeRaw.template(itemp,ichan,:))+(ichan)*maxheight,'color',cmap(itemp,:));
                end
            end
            grid on
            axis tight
            ax = axis;
            axis([ax(1), ax(2), 0, ax(4)]);
            title('Template');
            set(gca, 'YTickLabel', '');
            
            % write to PDF
            fig.Renderer='Painters';
            set(fig,'PaperOrientation','landscape');
            set(fig,'PaperUnits','normalized');
            set(fig,'PaperPosition', [0 0 1 1]);
            print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.label{ilabel},'_raster_all_templates_fromall.pdf']));
        end
        
    else
        fprintf('Could not find: %s\n',fname_spikes);
    end
    save(fname,'SpikeRaw','SpikeTrials');
end




