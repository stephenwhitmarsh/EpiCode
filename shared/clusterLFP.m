function [clusterindx, LFP_cluster] = clusterLFP(cfg, MuseStruct, force_clustering)

% Neccecary fields:
% cfg.cluster.name
cfg.visible             = ft_getopt(cfg, 'visible', 'on');
cfg.cluster.dbscan      = ft_getopt(cfg.cluster, 'dbscan', 'no');
cfg.cluster.kmeans      = ft_getopt(cfg.cluster, 'kmeans', 'no');
cfg.cluster.kmedoids    = ft_getopt(cfg.cluster, 'kmedoids', 'no');
if ~any([strcmp(cfg.cluster.dbscan, 'yes'), strcmp(cfg.cluster.kmeans, 'yes'), strcmp(cfg.cluster.kmedoids, 'yes')])
    fprintf('No cluster methods requested');
    return
end

cfg.cluster.reref           = ft_getopt(cfg.cluster, 'reref', 'no');
cfg.cluster.refmethod       = ft_getopt(cfg.cluster, 'refmethod', 'none');
cfg.cluster.ztransform      = ft_getopt(cfg.cluster, 'ztransform', 'no');
cfg.cluster.dbscan          = ft_getopt(cfg.cluster, 'dbscan', 'no');
cfg.cluster.epsilon         = ft_getopt(cfg.cluster, 'epsilon', 1500);
cfg.cluster.MinPts          = ft_getopt(cfg.cluster, 'MinPts', 10);
cfg.cluster.resamplefs      = ft_getopt(cfg.cluster, 'resamplefs', []);
cfg.cluster.channel         = ft_getopt(cfg.cluster, 'channel', 'all');
cfg.cluster.align           = ft_getopt(cfg.cluster, 'align', []);
cfg.cluster.align.latency   = ft_getopt(cfg.cluster.align, 'latency', 'all');

% load data if requested, and return
fname = fullfile(cfg.datasavedir,[cfg.prefix, 'clusterindx.mat']);

if nargin == 1
    if exist(fname, 'file')
        fprintf('Loading results clustering: %s\n', fname);
        % repeat to deal with load errors
        count = 0;
        err_count = 0;
        while count == err_count
            try
                load(fname, 'clusterindx', 'LFP_cluster');
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

if exist(fname, 'file') && force_clustering == false
    fprintf('Loading results clustering\n');
    load(fname, 'clusterindx', 'LFP_cluster');
    return
else
    clusterindx = [];
    LFP_cluster = [];
end

% read LFP according to LFP settings in configuration
cfg.LFP.write           = true;
cfg.LFP.name            = cfg.cluster.name;
LFP                     = readLFP(cfg, MuseStruct, force_clustering);

% preprocessing
for ipart = 1 : size(LFP, 2)

    for markername = string(fieldnames(LFP{ipart}))'

        if isempty(LFP{ipart}.(markername))
            continue
        end

        % select only those channels and time periods requested
        cfgtemp                 = [];
        cfgtemp.channel         = cfg.cluster.channel;
        LFP{ipart}.(markername) = ft_selectdata(cfgtemp, LFP{ipart}.(markername));

        % resample
        if ~isempty(cfg.cluster.resamplefs)
            cfgtemp                 = [];
            cfgtemp.resamplefs      = cfg.cluster.resamplefs;
            LFP{ipart}.(markername) = ft_resampledata(cfgtemp,LFP{ipart}.(markername));
        end

        % filtering and baseline correction
        cfgtemp                 = [];
        cfgtemp.lpfilter        = ft_getopt(cfg.cluster, 'lpfilter', 'no');
        cfgtemp.hpfilter        = ft_getopt(cfg.cluster, 'hpfilter', 'no');
        cfgtemp.bpfilter        = ft_getopt(cfg.cluster, 'bpfilter', 'no');
        cfgtemp.bsfilter        = ft_getopt(cfg.cluster, 'bsfilter', 'no');
        cfgtemp.lpfreq          = ft_getopt(cfg.cluster, 'lpfreq', []);
        cfgtemp.hpfreq          = ft_getopt(cfg.cluster, 'hpfreq', []);
        cfgtemp.bpfreq          = ft_getopt(cfg.cluster, 'bpfreq', []);
        cfgtemp.bsfreq          = ft_getopt(cfg.cluster, 'bsfreq', []);
        cfgtemp.demean          = ft_getopt(cfg.cluster, 'demean', 'no');
        cfgtemp.baselinewindow  = ft_getopt(cfg.cluster, 'baselinewindow', 'all');
        LFP{ipart}.(markername) = ft_preprocessing(cfgtemp, LFP{ipart}.(markername));

        % select only those time periods requested
        cfgtemp                 = [];
        cfgtemp.latency         = ft_getopt(cfg.cluster, 'latency', 'all');
        LFP{ipart}.(markername) = ft_selectdata(cfgtemp, LFP{ipart}.(markername));

        % bipolar rereferencing if requested
        if strcmp(cfg.cluster.reref, 'yes')
            if strcmp(cfg.cluster.refmethod, 'bipolar')
                labels_nonum    = regexprep(LFP{ipart}.(markername).label, '[0-9_]', '');
                [~,~,indx]      = unique(labels_nonum);
                for i = 1 : max(indx)
                    cfgtemp             = [];
                    cfgtemp.reref       = 'yes';
                    cfgtemp.refmethod   = 'bipolar';
                    cfgtemp.channel     = LFP{ipart}.(markername).label(indx==i);
                    group{i}            = ft_preprocessing(cfgtemp,LFP{ipart}.(markername));
                end
                LFP{ipart}.(markername)     = ft_appenddata([],group{:});
                clear group
            end
        end
    end
end

%%%%%%%%%%%%%%
%%% dbscan %%%
%%%%%%%%%%%%%%

if strcmp(cfg.cluster.dbscan, 'yes')
    for ipart = 1 : size(LFP, 2)

        for markername = string(fieldnames(LFP_avg))'

            if isempty(LFP{ipart}.(markername))
                continue
            end

            % put all in matrix with concatinated channels
            d = size(LFP{ipart}.trial{1});
            LFP_concatinated = nan(size(LFP{ipart}.(markername).trial,2), d(1)*d(2));
            for itrial = 1 : size(LFP{ipart}.(markername).trial, 2)
                LFP_concatinated(itrial,:) = reshape(LFP{ipart}.(markername).trial{itrial}', 1, d(1)*d(2));
            end

            % z-normalize data
            if strcmp(cfg.cluster.ztransform, 'yes')
                LFP_concatinated = nanznorm(LFP_concatinated);
                LFP_concatinated(isnan(LFP_concatinated)) = 0;
            end

            for k = 2 : 40
                [idx{k}, dist{k}] = knnsearch(LFP_concatinated, LFP_concatinated, 'K', k);
            end

            fig = figure('visible', cfg.visible); hold;
            for k = 2 : 40
                Y = sort(mean(dist{k}(:,2:end), 2));
                plot(Y);
            end
            grid on

            set(fig,'PaperOrientation','landscape');
            set(fig,'PaperUnits','normalized');
            set(fig,'PaperPosition', [0 0 1 1]);
            
            fname_fig = fullfile(cfg.imagesavedir, 'cluster', strcat(cfg.prefix, 'p', num2str(ipart), '_knn_', cfg.name{imarker}, '.png'));
            isdir_or_mkdir(fileparts(fname_fig));
            exportgraphics(fig, fname_fig);       
            
            close all

            [indx_dbscan, ~] = DBSCAN(LFP_concatinated, cfg.cluster.epsilon, cfg.cluster.MinPts);
            indx_dbscan = indx_dbscan + 1;

            clusterindx{ipart}.(markername).dbscan = indx_dbscan;

            fig = figure('visible', 'on'); hold;
            fig.Renderer = 'Painters';

            N = max(indx_dbscan);
            i = 1;
            for icluster = 1 : N

                % realign data again per cluster
                indx = indx_dbscan == icluster;

                [LFP_concatinated_realigned, nshift] = alignXcorr(LFP_concatinated(indx,:), 10);

                cfgtemp         = [];
                cfgtemp.trials  = find(indx)';
                LFP_shifted     = ft_selectdata(cfgtemp, LFP{ipart}.(markername));
                avg_original    = ft_timelockanalysis(cfgtemp, LFP{ipart}.(markername));

                for itrial = 1 : size(LFP_shifted.trial,2)
                    LFP_shifted.time{itrial} = LFP_shifted.time{itrial} + nshift(itrial) / LFP{ipart}.(markername).fsample;
                end

                cfgtemp         = [];
                avg_shifted     = ft_timelockanalysis(cfgtemp,LFP_shifted);

                subplot(3,N,icluster); hold;
                title(['Cluster ', num2str(icluster)]);
                imagesc(LFP_concatinated(indx_dbscan==icluster, :));
                set(gca,'xtick', []);
                set(gca,'fontsize', 6)
                axis tight

                subplot(3,N,icluster+N*1); hold;
                imagesc(LFP_concatinated_realigned);
                set(gca,'xtick', []);
                set(gca,'fontsize', 6)
                axis tight

                subplot(3,N,icluster+N*2); hold;
                h = max(max(abs(avg_shifted.avg)));
                n = 1;
                ytick = [];
                for ichan = 1 : size(avg_shifted.label,1)
                    ytick = [ytick, n*h];
                    plot(avg_original.time,avg_original.avg(ichan,:) + n*h, 'r');
                    plot(avg_shifted.time,avg_shifted.avg(ichan,:) + n*h, 'k');
                    temp = strsplit(avg_shifted.label{ichan}, '-');
                    avg_shifted.label{ichan} = temp{1};
                    n = n + 1;
                end
                yticks(ytick);
                set(gca,'TickLabelInterpreter', 'none')
                set(gca,'fontsize', 6)
                yticklabels(avg_shifted.label);
                xlabel('Time (s)');
                axis tight
                xlim([-0.5 1]);
                box off

                LFP_cluster{ipart}.(markername).dbscan{icluster} = sel;
            end

            
            set(fig,'PaperOrientation','landscape');
            set(fig,'PaperUnits','normalized');
            set(fig,'PaperPosition', [0 0 1 1]);
            fname_fig = fullfile(cfg.imagesavedir, 'cluster', strcat(cfg.prefix, 'p', num2str(ipart), '_DBSCAN_', markername, '_LFP.png'));
            isdir_or_mkdir(fileparts(fname_fig));
            exportgraphics(fig, fname_fig);
            close all

        end
    end
end

%%%%%%%%%%%%%%
%%% kmeans %%%
%%%%%%%%%%%%%%

if strcmp(cfg.cluster.kmeans, 'yes')
    for ipart = 1 : size(LFP,2)
        for imarker = 1 : size(LFP{ipart},2)

            if isempty(LFP{ipart}.(markername))
                continue
            end            % put all in matrix with concatinated channels
            d = size(LFP{ipart}.(markername).trial{1});
            LFP_concatinated = nan(size(LFP{ipart}.(markername).trial,2),d(1)*d(2));
            for itrial = 1 : size(LFP{ipart}.(markername).trial,2)
                LFP_concatinated(itrial,:) = reshape(LFP{ipart}.(markername).trial{itrial}',1,d(1)*d(2));
            end

            % z-normalize data
            if strcmp(cfg.cluster.ztransform, 'yes')
                LFP_concatinated = nanznorm(LFP_concatinated);
                LFP_concatinated(isnan(LFP_concatinated)) = 0;
            end

            eva = evalclusters(LFP_concatinated,'kmeans','CalinskiHarabasz','KList',[1:10]);
            indx_kmeans = eva.OptimalY;
%             for N = cfg.cluster.N

%                 [indx_kmeans, ~] = kmeans(LFP_concatinated, N, 'MaxIter', 2000, 'Replicates', 100, 'dist', 'sqeuclidean');

                clusterindx{ipart}.(markername).kmeans = indx_kmeans;

                % plot kmeans clusters
                fig = figure('visible', cfg.visible);
                fig.Renderer = 'Painters';

                for icluster = 1 : eva.OptimalK

                    % realign data again per cluster
                    indx = indx_kmeans == icluster;

                    [LFP_concatinated_realigned, nshift] = alignXcorr(LFP_concatinated(indx,:), 10);

                    cfgtemp         = [];
                    cfgtemp.trials  = find(indx)';
                    LFP_shifted     = ft_selectdata(cfgtemp, LFP{ipart}.(markername));
                    avg_original    = ft_timelockanalysis(cfgtemp, LFP{ipart}.(markername));

                    for itrial = 1 : size(LFP_shifted.trial,2)
                        LFP_shifted.time{itrial} = LFP_shifted.time{itrial} + nshift(itrial) / LFP{ipart}.(markername).fsample;
                    end

                    cfgtemp         = [];
                    avg_shifted     = ft_timelockanalysis(cfgtemp,LFP_shifted);

                    subplot(3, eva.OptimalK, icluster); hold;
                    title(['Cluster ', num2str(icluster-1)]);

                    imagesc(LFP_concatinated(indx_kmeans==icluster, :));
                    set(gca,'xtick', []);
                    set(gca,'fontsize', 6)
                    axis tight

                    subplot(3, eva.OptimalK, icluster+eva.OptimalK*1); hold;
                    imagesc(LFP_concatinated_realigned);
                    set(gca,'xtick', []);
                    set(gca,'fontsize', 6)
                    axis tight

                    subplot(3,eva.OptimalK,icluster+eva.OptimalK*2); hold;
                    h = max(max(abs(avg_shifted.avg)));
                    n = 1;
                    ytick = [];
                    for ichan = 1 : size(avg_shifted.label,1)
                        ytick = [ytick, n*h];
                        plot(avg_original.time,avg_original.avg(ichan,:) + n*h, 'r');
                        plot(avg_shifted.time,avg_shifted.avg(ichan,:) + n*h, 'k');
                        temp = strsplit(avg_shifted.label{ichan}, '-');
                        label{ichan} = temp{1};
                        n = n + 1;
                    end
                    yticks(ytick);
                    set(gca,'TickLabelInterpreter', 'none')
                    set(gca,'fontsize', 6)
                    yticklabels(label);
                    xlabel('Time (s)');
                    axis tight
                    xlim([LFP{ipart}.(markername).time{1}(1), LFP{ipart}.(markername).time{1}(end)]);
                    ax = axis;
                    if ~isstring(cfg.cluster.align.latency)
                        patch([cfg.cluster.align.latency(1), cfg.cluster.align.latency(2), ...
                            cfg.cluster.align.latency(2), cfg.cluster.align.latency(1)], ...
                            [ax(3), ax(3), ax(4), ax(4)], 'r', 'facealpha', 0.1, 'edgecolor', 'none');
                    end
                    box off

                    % save for output
                    LFP_cluster{ipart}.(markername).kmeans{icluster} = avg_shifted;
                end

                set(fig,'PaperOrientation', 'landscape');
                set(fig,'PaperUnits', 'normalized');
                set(fig,'PaperPosition', [0 0 1 1]);
                fname_fig = fullfile(cfg.imagesavedir, 'cluster', strcat(cfg.prefix, 'p' ,num2str(ipart), '_kmeans_', markername, '_N', num2str(eva.OptimalK), '_LFP.png'));
                isdir_or_mkdir(fileparts(fname_fig));
                exportgraphics(fig, fname_fig);
                
                close all
                %
                %                 % plot silhouette
                %                 fig = figure('visible', 'off'); hold;
                %                 [silh,h] = silhouette(LFP_concatinated, indx_kmeans, 'sqeuclidean');
                %                 set(fig,'PaperOrientation', 'landscape');
%                 set(fig,'PaperUnits', 'normalized');
%                 set(fig,'PaperPosition', [0 0 1 1]);
%                 print(fig, '-dpdf', fullfile(cfg.imagesavedir, [cfg.prefix, 'p' ,num2str(ipart), '_kmeans_', cfg.name{imarker}, '_N', num2str(N), '_silhouette.pdf']));
%                 print(fig, '-dpng', fullfile(cfg.imagesavedir, [cfg.prefix, 'p' ,num2str(ipart), '_kmeans_', cfg.name{imarker}, '_N', num2str(N), '_silhouette.png']));
%                 close all
%             end
        end
    end
end

%%%%%%%%%%%%%%%%
%%% kmedoids %%%
%%%%%%%%%%%%%%%%

if strcmp(cfg.cluster.kmedoids, 'yes')
    for ipart = 1 : size(LFP,2)
        for imarker = 1 : size(LFP{ipart},2)

            if isempty(LFP{ipart}.(markername))
                continue
            end

            % put all in matrix with concatinated channels
            d = size(LFP{ipart}.(markername).trial{1});
            LFP_concatinated = nan(size(LFP{ipart}.(markername).trial,2),d(1)*d(2));
            for itrial = 1 : size(LFP{ipart}.(markername).trial,2)
                LFP_concatinated(itrial,:) = reshape(LFP{ipart}.(markername).trial{itrial}',1,d(1)*d(2));
            end

            % z-normalize data
            if strcmp(cfg.cluster.ztransform, 'yes')
                LFP_concatinated = nanznorm(LFP_concatinated);
                LFP_concatinated(isnan(LFP_concatinated)) = 0;
            end

            for N = cfg.cluster.N

                [indx_kmedoids, ~] = kmedoids(LFP_concatinated, N);
                clusterindx{ipart}.(markername).kmedoids(:,N) = indx_kmedoids;

                % plot kmedoids clusters
                fig = figure('visible', cfg.visible);
                fig.Renderer = 'Painters';

                for icluster = 1 : N

                    % realign data again per cluster
                    indx = indx_kmedoids == icluster;

                    % select period for alignemnt
                    cfgtemp         = [];
                    cfgtemp.latency = cfg.cluster.align.latency;
                    cfgtemp.trials  = indx;
                    LFP_temp        = ft_selectdata(cfgtemp, LFP{ipart}.(markername));

                    % put all in matrix with concatinated channels
                    d = size(LFP_temp.trial{1});
                    LFP_sel_concatinated = nan(size(LFP_temp.trial,2),d(1)*d(2));
                    for itrial = 1 : size(LFP_temp.trial,2)
                        LFP_sel_concatinated(itrial,:) = reshape(LFP_temp.trial{itrial}', 1, d(1)*d(2));
                    end

                    % align
                    [LFP_concatinated_realigned, nshift] = alignXcorr(LFP_sel_concatinated, 10);

                    % create averages
                    cfgtemp         = [];
                    cfgtemp.trials  = find(indx)';
                    LFP_shifted     = ft_selectdata(cfgtemp, LFP{ipart}.(markername));
                    avg_original    = ft_timelockanalysis(cfgtemp, LFP{ipart}.(markername));

                    for itrial = 1 : size(LFP_shifted.trial,2)
                        LFP_shifted.time{itrial} = LFP_shifted.time{itrial} + nshift(itrial) / LFP{ipart}.(markername).fsample;
                    end

                    cfgtemp         = [];
                    avg_shifted     = ft_timelockanalysis(cfgtemp,LFP_shifted);

                    % plot matrix pre-alignment
                    subplot(3, N, icluster); hold;
                    title(['Cluster ', num2str(icluster)]);
                    imagesc(LFP_sel_concatinated);
                    set(gca,'xtick', []);
                    set(gca,'fontsize', 6)
                    axis tight

                    % plot matrix post-alignment
                    subplot(3, N, icluster+N*1); hold;
                    imagesc(LFP_concatinated_realigned);
                    set(gca,'xtick', []);
                    set(gca,'fontsize', 6)
                    axis tight

                    % plot LFP pre- and post-alignment
                    subplot(3, N, icluster+N*2); hold;
                    h = max(max(abs(avg_shifted.avg)));
                    n = 1;
                    ytick = [];
                    for ichan = 1 : size(avg_shifted.label,1)
                        ytick = [ytick, n*h];
                        plot(avg_original.time, avg_original.avg(ichan,:) + n*h, 'r');
                        plot(avg_shifted.time, avg_shifted.avg(ichan,:) + n*h, 'k');
                        temp = strsplit(avg_shifted.label{ichan}, '-');
                        label{ichan} = temp{1};
                        n = n + 1;
                    end
                    yticks(ytick);
                    set(gca,'TickLabelInterpreter', 'none')
                    set(gca,'fontsize', 6)
                    yticklabels(label);
                    xlabel('Time (s)');
                    axis tight
                    xlim([LFP{ipart}.(markername).time{1}(1), LFP{ipart}.(markername).time{1}(end)]);
                    ax = axis;
                    if ~isstring(cfg.cluster.align.latency)
                        patch([cfg.cluster.align.latency(1), cfg.cluster.align.latency(2), ...
                            cfg.cluster.align.latency(2), cfg.cluster.align.latency(1)], ...
                            [ax(3), ax(3), ax(4), ax(4)], 'r', 'facealpha', 0.1, 'edgecolor', 'none');
                    end
                    box off

                    % save for output
                    LFP_cluster{ipart}.(markername).kmedoids{N}{icluster} = avg_shifted;
                end

                set(fig,'PaperOrientation', 'landscape');
                set(fig,'PaperUnits', 'normalized');
                set(fig,'PaperPosition', [0 0 1 1]);            
                fname_fig = fullfile(cfg.imagesavedir, 'cluster', strcat(cfg.prefix, 'p' ,num2str(ipart), '_kmedoids_', markername, '_N', num2str(N), '_LFP.png'));
                isdir_or_mkdir(fileparts(fname_fig));
                exportgraphics(fig, fname_fig);
                
                close all

%                 % plot silhouette
%                 fig = figure('visible', 'off'); hold;
%                 [~, ~] = silhouette(LFP_concatinated, indx_kmedoids, 'sqeuclidean');
%                 set(fig,'PaperOrientation', 'landscape');
%                 set(fig,'PaperUnits', 'normalized');
%                 set(fig,'PaperPosition', [0 0 1 1]);
%                 print(fig, '-dpdf', fullfile(cfg.imagesavedir, strcat(cfg.prefix, 'p' ,num2str(ipart), '_kmedoids_', markername, '_N', num2str(N)), '_silhouette.pdf'));
%                 print(fig, '-dpng', fullfile(cfg.imagesavedir, strcat(cfg.prefix, 'p' ,num2str(ipart), '_kmedoids_', cfg.name{imarker}, '_N', num2str(N)), '_silhouette.png'));
%                 close all
            end
        end
    end
end

save(fname, 'clusterindx', 'LFP_cluster');


%                 %% Dynamic Time Warping
%
%                 maxamp = 100;
%                 tic;
%
%                 if maxamp == 0
%                     dtwD = pdist(z_input,@(Xi,Xj) dtwdist(Xi,Xj,'squared'));
%                 else
%                     dtwD = pdist(z_input,@(Xi,Xj) dtwdist(Xi,Xj,maxamp,'squared'));
%                 end
%
%                 save(fullfile(cfg.datasavedir,[cfg.prefix,'dtwD_sel.mat']),'dtwD');
%                 elapsedTime = toc
%
%                 fig = figure;
%                 fig.Renderer = 'Painters';
%                 imagesc(squareform(dtwD));
%                 colorbar;
%                 title('Distance Matrix');
%                 set(fig,'PaperOrientation','landscape');
%                 set(fig,'PaperUnits','normalized');
%                 set(fig,'PaperPosition', [0 0 1 1]);
%                 print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'distance_matrix_DTW_',cfg.name{imarker},'_N',num2str(N),'.pdf']),'-r600');
%                 print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,'distance_matrix_DTW_',cfg.name{imarker},'_N',num2str(N),'.png']),'-r600');
%                 close all
%
%
%                 % plot tree
%
%                 clustTreeEuc = linkage(dtwD,'weighted');
%                 clear sumc
%
%                 fig = figure;
%                 fig.Renderer = 'Painters';
%
%                 subplot(1,3,1);
%                 [h,nodes] = dendrogram(clustTreeEuc,0);
%                 title(sprintf('Cophenetic corr: %0.2f',cophenet(clustTreeEuc,eucD)));
%                 set(gca,'xtick',[]);
%
%                 subplot(1,3,2);
%                 [h,nodes] = dendrogram(clustTreeEuc,10);
%                 title('Max leafs: 10');
%
%                 subplot(1,3,3);
%                 [h,nodes] = dendrogram(clustTreeEuc,N);
%                 title(sprintf('Max leafs: %d',N));
%                 xlabel('Leaf number');
%
%                 for i = 1 : N
%                     sumc(i) = sum(nodes==i);
%                 end
%                 set(gca,'Xticklabel',sumc);
%                 xlabel('#Datapoints');
%
%                 %         for i = 1 : N
%                 %             cluster{i}          = input(nodes==i,:);
%                 %             cluster_mean(i,:)   = mean(cluster{i});
%                 %         end
%                 %         ymin = min(min(cluster_mean));
%                 %         ymax = max(max(cluster_mean));
%                 %
%
%                 set(fig,'PaperOrientation','landscape');
%                 set(fig,'PaperUnits','normalized');
%                 set(fig,'PaperPosition', [0 0 1 1]);
%                 print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'dendogram_DTW_',cfg.name{imarker},'.pdf']),'-r600');
%                 print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,'dendogram_DTW_',cfg.name{imarker},'.png']),'-r600');
%                 close all
%
%                 % plot grids
%                 fig = figure;
%                 fig.Renderer = 'Painters';
%
%                 for i = 1 : N
%
%                     [z_input_realigned,~] = align_data(z_input(indx==i,:),5);
%
%                     subplot(3,N,i);
%                     imagesc(z_input(indx==i,:));
%                     set(gca,'xtick',[]);
%
%                     subplot(3,N,i+N*1);
%                     imagesc(z_input_realigned);
%                     set(gca,'xtick',[]);
%
%                     subplot(3,N,i+N*2); hold;
%                     plot(dat{imarker}.time{1},nanmean(z_input(indx==i,:)));
%                     %             xlim([dat{imarker}.time{1}(1),dat{imarker}.time{1}(end)]);
%                     %             ylim([ymin*1.1,ymax*1.1]);
%                     axis tight
%                     plot(dat{imarker}.time{1},nanmean(z_input_realigned));
%
%                 end
%                 %         subplot(3,N,N*2+1);
%                 %         [silh,h] = silhouette(shifted_clean_z,nodes,'euclidean');
%
%                 set(fig,'PaperOrientation','landscape');
%                 set(fig,'PaperUnits','normalized');
%                 set(fig,'PaperPosition', [0 0 1 1]);
%                 print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'hierach_clusters_DTW_',cfg.name{imarker},'_N',num2str(N),'.pdf']),'-r600');
%                 print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,'hierach_clusters_DTW_',cfg.name{imarker},'_N',num2str(N),'.png']),'-r600');
%                 close all
%
