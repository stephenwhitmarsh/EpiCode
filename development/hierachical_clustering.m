
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
            