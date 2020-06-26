function [clusterindx] = findPattern2D(cfg, MuseStruct, force_clustering)

fname = fullfile(cfg.datasavedir,[cfg.prefix,'clusterindx.mat']);
if exist(fname,'file') && force_clustering == false
    fprintf('*******************************\n');
    fprintf('** Loading results alignment **\n');
    fprintf('*******************************\n\n');
    load(fname,'clusterindx');
    return
end

% read LFP according to LFP settings in configuration
cfgtemp             = rmfield(cfg,'LFP');
cfgtemp.LFP.name    = cfg.cluster.name;
cfgtemp.LFP.channel = cfg.cluster.channel;
cfgtemp.LFP.write   = false;
[LFP]               = readLFP(cfgtemp, MuseStruct, true);

% bipolar rereferencing if requested
cfg.cluster.reref       = ft_getopt(cfg.cluster,'reref','no');
cfg.cluster.refmethod   = ft_getopt(cfg.cluster,'refmethod','none');
if strcmp(cfg.cluster.refmethod, 'bipolar')
    for ipart = 1:size(MuseStruct,2)
        for imarker = 1 : size(cfg.cluster.name,1)
            labels_nonum    = regexprep(LFP{ipart}{imarker}.label, '[0-9_]', '');
            [~,~,indx]      = unique(labels_nonum);
            clear group
            for i = 1 : max(indx)
                cfgtemp             = [];
                cfgtemp.reref       = 'yes';
                cfgtemp.refmethod   = 'bipolar';
                cfgtemp.channel     = LFP{ipart}{imarker}.label(indx==i);
                group{i}            = ft_preprocessing(cfgtemp,LFP{ipart}{imarker});
            end
            LFP{ipart}{imarker} = ft_appenddata([],group{:});
            clear group
        end
    end
end

% cluster data per part
for ipart = 1 : size(cfg.directorylist,2)
    for imarker = 1 : size(cfg.cluster.name,1)
        
        % baseline correction
        cfgtemp                 = [];
        cfgtemp.demean          = 'yes';
        cfgtemp.baselinewindow  = [-0.5 -0.2];
        LFP{ipart}{imarker}     = ft_preprocessing(cfgtemp,LFP{ipart}{imarker});
        
        % baseline correction
        cfgtemp                 = [];
        cfgtemp.latency         = [-0.5 1];
        LFP{ipart}{imarker}     = ft_selectdata(cfgtemp,LFP{ipart}{imarker});
        
        % put all in matrix (concatinated)
        d = size(LFP{ipart}{imarker}.trial{1});
        LFP_concatinated = nan(size(LFP{ipart}{imarker}.trial,2),d(1)*d(2));
        for itrial = 1 : size(LFP{ipart}{imarker}.trial,2)
            LFP_concatinated(itrial,:) = reshape(LFP{ipart}{imarker}.trial{itrial}',1,d(1)*d(2));
        end
        clear temp
        
        % z-normalize data
        LFP_concatinated_z = nanznorm(LFP_concatinated);
        LFP_concatinated_z(isnan(LFP_concatinated_z)) = 0;
        
        %%%%%%%%%%%%%%
        %%% dbscan %%%
        %%%%%%%%%%%%%%
        
        for k = 2 : 40
            [idx{k}, dist{k}] = knnsearch(LFP_concatinated,LFP_concatinated,'K',k);
        end
        
        fig = figure('visible','off'); hold;
        for k = 2 : 40
            Y = sort(mean(dist{k}(:,2:end),2));
            plot(Y);
        end
        grid on
        
        set(fig,'PaperOrientation','landscape');
        set(fig,'PaperUnits','normalized');
        set(fig,'PaperPosition', [0 0 1 1]);
        print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'p',num2str(ipart),'_knn_',cfg.name{imarker},'.pdf']),'-r600');
        print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,'p',num2str(ipart),'_knn_',cfg.name{imarker},'.png']),'-r600');
        close all
        
        MinPts  = 10;
        epsilon = 4000;
        [indx_dbscan{ipart}{imarker}, ~] = DBSCAN(LFP_concatinated,epsilon,MinPts);
        indx_dbscan{ipart}{imarker} = indx_dbscan{ipart}{imarker} + 1;
        
        fig = figure('visible','on'); hold;
        fig.Renderer = 'Painters';
        
        N = max(indx_dbscan{ipart}{imarker});
        i = 1;
        for icluster = 1 : N
            
            indx = indx_dbscan{ipart}{imarker} == icluster;
            [LFP_concatinated_realigned, nshift] = alignXcorr(LFP_concatinated(indx,:),10);
            
            LFP_shifted = LFP{ipart}{imarker};
            for itrial = find(indx_dbscan{ipart}{imarker}==icluster)'
                LFP_shifted.time{itrial} = LFP{ipart}{imarker}.time{itrial} + nshift(indx(itrial)) / LFP{ipart}{imarker}.fsample;
            end
            cfgtemp         = [];
            cfgtemp.trials  = find(indx_dbscan{ipart}{imarker}==icluster);
            sel             = ft_timelockanalysis(cfgtemp,LFP_shifted);
            
            subplot(3,N,i); hold;
            if icluster > 0
                title(['Cluster ', num2str(i)]);
            else
                title('Noise');
            end
            
            imagesc(LFP_concatinated(indx_dbscan{ipart}{imarker}==icluster,:));
            set(gca,'xtick',[]);
            set(gca,'fontsize',6)
            axis tight
            
            subplot(3,N,i+N*1); hold;
            imagesc(LFP_concatinated_realigned);
            set(gca,'xtick',[]);
            set(gca,'fontsize',6)
            axis tight
            
            subplot(3,N,i+N*2); hold;
            h = max(max(abs(sel.avg)));
            n = 1;
            ytick = [];
            for ichan = 1 : size(sel.label,1)
                ytick = [ytick, n*h];
                plot(sel.time,sel.avg(ichan,:) + n*h,'color','k');
                temp = strsplit(sel.label{ichan},'-');
                sel.label{ichan} = temp{1};
                n = n + 1;
            end
            yticks(ytick);
            set(gca,'TickLabelInterpreter','none')
            set(gca,'fontsize',6)
            yticklabels(sel.label);
            xlabel('Time (s)');
            axis tight
            xlim([-0.5 1]);
            box off
            
            i = i + 1;
        end
        
        set(fig,'PaperOrientation','landscape');
        set(fig,'PaperUnits','normalized');
        set(fig,'PaperPosition', [0 0 1 1]);
        print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'p',num2str(ipart),'_kmediods_',cfg.name{imarker},'.pdf']),'-r600');
        print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,'p',num2str(ipart),'_kmediods_',cfg.name{imarker},'.png']),'-r600');
        close all
        
        %% Euclidian distance for hierarchical clustering
        %             eucD = pdist(LFP_concatinated_z,'seuclidean');
        %             eucD = pdist(shifted_clean_z,@naneucdist); % method to allow nans
        %
        %             [Y,stress,disparities] = mdscale(eucD,3,'Criterion','sstress');
        %             figure;
        %             scatter3(Y(:,1),Y(:,2),Y(:,3));
        
        %% Hierarchical Clustering
        %             clustTreeEuc = linkage(eucD,'average');
        
        % loop over increasing nr of clusters
        for N = 2 : 6
            
            %                 % plot tree Hierarchical Clustering
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
            %                 for i = 1 : N
            %                     cluster{i}          = LFP_concatinated(nodes==i,:);
            %                     cluster_mean(i,:)   = mean(cluster{i});
            %                 end
            %                 ymin = min(min(cluster_mean));
            %                 ymax = max(max(cluster_mean));
            %
            %                 set(fig,'PaperOrientation','landscape');
            %                 set(fig,'PaperUnits','normalized');
            %                 set(fig,'PaperPosition', [0 0 1 1]);
            %                 %                 print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'dendogram_',cfg.name{imarker},'_N',num2str(N),'.pdf']),'-r600');
            %                 print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,'dendogram_',cfg.name{imarker},'_N',num2str(N),'.png']));
            %                 close all
            %
            %
            %                 % plot data
            %
            %                 fig = figure;
            %                 fig.Renderer = 'Painters';
            %
            %                 for i = 1 : N
            %
            %                     if size(find(nodes==i),1) > 1
            %                         [LFP_concatinated_aligned,~] = align_data(LFP_concatinated_z(nodes==i,:),5);
            %                     else
            %                         LFP_concatinated_aligned = LFP_concatinated_z(nodes==i,:);
            %                     end
            %
            %                     subplot(3,N,i);
            %                     imagesc(LFP_concatinated_z(nodes==i,:));
            %                     set(gca,'xtick',[]);
            %
            %                     subplot(3,N,i+N*1);
            %                     imagesc(LFP_concatinated_aligned);
            %                     set(gca,'xtick',[]);
            %
            %                     subplot(3,N,i+N*2); hold;
            %
            %                     cfgtemp         = [];
            %                     cfgtemp.trials  = find(nodes == i);
            %                     sel             = ft_timelockanalysis(cfgtemp,LFP{ipart}{imarker});
            %
            %                     plot(sel.time, sel.avg,'color',[0.7, 0.7, 0.7]);
            %                     plot(sel.time, mean(sel.avg),'color','k');
            %
            %                 end
            %                 %         subplot(3,N,N*2+1);
            %                 %         [silh,h] = silhouette(shifted_clean_z,nodes,'euclidean');
            %                 %
            %                 set(fig,'PaperOrientation','landscape');
            %                 set(fig,'PaperUnits','normalized');
            %                 set(fig,'PaperPosition', [0 0 1 1]);
            %                 % print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'hierach_clusters_',cfg.name{imarker},'_N',num2str(N),'.pdf']),'-r600');
            %                 print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,'hierach_clusters_',cfg.name{imarker},'_N',num2str(N),'.png']));
            %                 close all
            %                 dbscan
            
            %%%%%%%%%%%%%%
            %%% kmeans %%%
            %%%%%%%%%%%%%%
            
            [indx_kmeans{ipart}{imarker}(:,N), ~] = kmeans(LFP_concatinated,N,'MaxIter',2000,'Replicates',100,'dist','sqeuclidean');
            
            clear cluster cluster_mean
            for i = 1 : N
                cluster{i}          = LFP_concatinated(indx_kmeans{ipart}{imarker}(:,N)==i,:);
                cluster_mean(i,:)   = mean(cluster{i});
            end
            ymin = min(min(cluster_mean))*1.3;
            ymax = max(max(cluster_mean))*1.3;
            
            % plot kmeans clusters
            fig = figure(1); hold;
            fig.Renderer = 'Painters';
            
            for i = 1 : N
                
                % realign data again, but not per per cluster
                indx = indx_kmeans{ipart}{imarker}(:,N) == i;
                [LFP_concatinated_realigned, nshift] = alignXcorr(LFP_concatinated(indx,:),10);
                
                LFP_shifted = LFP{ipart}{imarker};
                for itrial = find(indx_kmeans{ipart}{imarker}(:,N)==i)'
                    LFP_shifted.time{itrial} = LFP{ipart}{imarker}.time{itrial} + nshift(indx(itrial)) / LFP{ipart}{imarker}.fsample;
                end
                cfgtemp         = [];
                cfgtemp.trials  = find(indx_kmeans{ipart}{imarker}(:,N)==i);
                sel             = ft_timelockanalysis(cfgtemp,LFP_shifted);
                
                subplot(3,N,i); hold;
                title(['Cluster ', num2str(i)]);
                imagesc(LFP_concatinated(indx_kmeans{ipart}{imarker}(:,N)==i,:));
                set(gca,'xtick',[]);
                set(gca,'fontsize',6)
                axis tight
                
                subplot(3,N,i+N*1); hold;
                imagesc(LFP_concatinated_realigned);
                set(gca,'xtick',[]);
                set(gca,'fontsize',6)
                axis tight
                
                subplot(3,N,i+N*2); hold;
                h = max(max(abs(sel.avg)));
                n = 1;
                ytick = [];
                for ichan = 1 : size(sel.label,1)
                    ytick = [ytick, n*h];
                    plot(sel.time,sel.avg(ichan,:) + n*h,'color','k');
                    temp = strsplit(sel.label{ichan},'-');
                    sel.label{ichan} = temp{1};
                    n = n + 1;
                end
                yticks(ytick);
                set(gca,'TickLabelInterpreter','none')
                set(gca,'fontsize',6)
                yticklabels(sel.label);
                xlabel('Time (s)');
                axis tight
                xlim([-0.5 1]);
                box off
            end
            
            set(fig,'PaperOrientation','landscape');
            set(fig,'PaperUnits','normalized');
            set(fig,'PaperPosition', [0 0 1 1]);
            print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'p',num2str(ipart),'_kmeans_',cfg.name{imarker},'_N',num2str(N),'.pdf']),'-r600');
            print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,'p',num2str(ipart),'_kmeans_',cfg.name{imarker},'_N',num2str(N),'.png']),'-r600');
            close all
            
            % plot siluette
            fig = figure(2);
            [silh,h] = silhouette(LFP_concatinated,indx_kmeans{ipart}{imarker}(:,N),'sqeuclidean');
            set(fig,'PaperOrientation','landscape');
            set(fig,'PaperUnits','normalized');
            set(fig,'PaperPosition', [0 0 1 1]);
            print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'p',num2str(ipart),'_kmeans_silhouette_',cfg.name{imarker},'_N',num2str(N),'.pdf']));
            print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,'p',num2str(ipart),'_kmeans_silhouette_',cfg.name{imarker},'_N',num2str(N),'.png']));
            close all
            
            %%%%%%%%%%%%%%%%
            %%% kmediods %%%
            %%%%%%%%%%%%%%%%
            
            [indx_kmediods{ipart}{imarker}(:,N), ~] = kmedoids(LFP_concatinated,N);
            
            clear cluster cluster_mean
            for i = 1 : N
                cluster{i}          = LFP_concatinated(indx_kmediods{ipart}{imarker}(:,N)==i,:);
                cluster_mean(i,:)   = mean(cluster{i});
            end
            ymin = min(min(cluster_mean))*1.3;
            ymax = max(max(cluster_mean))*1.3;
            
            % plot kmediods clusters
            fig = figure('visible','off'); hold;
            fig.Renderer = 'Painters';
            
            for i = 1 : N
                
                % realign data again, but not per per cluster
                indx = indx_kmediods{ipart}{imarker}(:,N) == i;
                [LFP_concatinated_realigned, nshift] = alignXcorr(LFP_concatinated(indx,:),10);
                
                LFP_shifted = LFP{ipart}{imarker};
                for itrial = find(indx_kmediods{ipart}{imarker}(:,N)==i)'
                    LFP_shifted.time{itrial} = LFP{ipart}{imarker}.time{itrial} + nshift(indx(itrial)) / LFP{ipart}{imarker}.fsample;
                end
                cfgtemp         = [];
                cfgtemp.trials  = find(indx_kmediods{ipart}{imarker}(:,N)==i);
                sel             = ft_timelockanalysis(cfgtemp,LFP_shifted);
                
                subplot(3,N,i); hold;
                title(['Cluster ', num2str(i)]);
                imagesc(LFP_concatinated(indx_kmediods{ipart}{imarker}(:,N)==i,:));
                set(gca,'xtick',[]);
                set(gca,'fontsize',6)
                axis tight
                
                subplot(3,N,i+N*1); hold;
                imagesc(LFP_concatinated_realigned);
                set(gca,'xtick',[]);
                set(gca,'fontsize',6)
                axis tight
                
                subplot(3,N,i+N*2); hold;
                h = max(max(abs(sel.avg)));
                n = 1;
                ytick = [];
                for ichan = 1 : size(sel.label,1)
                    ytick = [ytick, n*h];
                    plot(sel.time,sel.avg(ichan,:) + n*h,'color','k');
                    temp = strsplit(sel.label{ichan},'-');
                    sel.label{ichan} = temp{1};
                    n = n + 1;
                end
                yticks(ytick);
                set(gca,'TickLabelInterpreter','none')
                set(gca,'fontsize',6)
                yticklabels(sel.label);
                xlabel('Time (s)');
                axis tight
                xlim([-0.5 1]);
                box off
            end
            
            set(fig,'PaperOrientation','landscape');
            set(fig,'PaperUnits','normalized');
            set(fig,'PaperPosition', [0 0 1 1]);
            print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'p',num2str(ipart),'_kmediods_',cfg.name{imarker},'_N',num2str(N),'.pdf']),'-r600');
            print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,'p',num2str(ipart),'_kmediods_',cfg.name{imarker},'_N',num2str(N),'.png']),'-r600');
            close all
            
            % plot siluette
            fig = figure('visible','off');
            [silh,h] = silhouette(LFP_concatinated,indx_kmediods{ipart}{imarker}(:,N),'sqeuclidean');
            set(fig,'PaperOrientation','landscape');
            set(fig,'PaperUnits','normalized');
            set(fig,'PaperPosition', [0 0 1 1]);
            print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'p',num2str(ipart),'_kmediods_silhouette_',cfg.name{imarker},'_N',num2str(N),'.pdf']));
            print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,'p',num2str(ipart),'_kmediods_silhouette_',cfg.name{imarker},'_N',num2str(N),'.png']));
            close all
            
            
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
        end
    end
end

clusterindx{ipart}{imarker}.kmeans      = indx_kmeans;
clusterindx{ipart}{imarker}.kmediods    = indx_kmediods;
clusterindx{ipart}{imarker}.dbscan      = indx_dbscan;
save(fname,'clusterindx');
