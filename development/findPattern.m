function findPattern(cfg, dat_in, force)

fname = fullfile(cfg.datasavedir,[cfg.prefix,'patterns_found.mat']);
if exist(fname,'file') && force == false
    fprintf('*******************************\n');
    fprintf('** Loading results alignment **\n');
    fprintf('*******************************\n\n');
    load(fullfile(cfg.datasavedir,[cfg.prefix,'patterns_found.mat']),'patterns_found');
else
    
    
    
    for imarker = 1 : size(dat_in,2)
        
        
        if force == true
            fprintf('*********************************\n');
            fprintf('** Forced redoing of alignment **\n');
            fprintf('*********************************\n\n');
        else
            fprintf('**************\n');
            fprintf('** Aligning **\n');
            fprintf('**************\n\n');
        end
        
        % baseline correct
        cfgtemp = [];
        cfgtemp.demean = 'yes';
        cfgtemp.baselinewindow = [-0.6 -0.3];
        dat{imarker} = ft_preprocessing(cfgtemp,dat_in{imarker});
        
        % select data
        cfgtemp = [];
        cfgtemp.channel = 1;
        cfgtemp.latency = [-0.1, 0.2]; % for alignement best is -0.1:0.2
        dat{imarker}  = ft_selectdata(cfgtemp,dat{imarker});
        
        % put all in matrix
        clear input
        for itrial = 1 : size(dat{imarker}.trial,2)
            input(itrial,:) = dat{imarker}.trial{itrial}(1,:);
        end
        
        % align data with closest-peak in cross correlation
        [shifted,nshift] = align_data(input,5);
        
        % remove those that moved too much
        rejected = (nshift < mean(nshift)-std(nshift)*2) | (nshift > mean(nshift) + std(nshift)*2);
        shifted_clean = shifted(~rejected,:);
        
        % z-normalize data
        shifted_clean_z = nanznorm(shifted_clean);
        shifted_clean_z(isnan(shifted_clean_z)) = 0;
        
        fig = figure;
        fig.Renderer = 'Painters';

        subplot(2,5,1);
        imagesc(input(~rejected,:));
        title('Cleaned');
        set(gca,'xtick',[]);
        subplot(2,5,6);
        plot(dat{imarker}.time{1},input(~rejected,:)');
        axis tight
        
        subplot(2,5,2);
        imagesc(input(rejected,:));
        title('Rejected');
        set(gca,'xtick',[]);
        subplot(2,5,7);
        plot(dat{imarker}.time{1},input(rejected,:)');
        axis tight
        
        subplot(2,5,3);
        imagesc(shifted(~rejected,:));
        title('Cleaned & Aligned');
        set(gca,'xtick',[]);
        subplot(2,5,8);
        plot(dat{imarker}.time{1},shifted(~rejected,:)');
        axis tight
        
        subplot(2,5,4);
        imagesc(shifted(rejected,:));
        title('Misaligned');
        set(gca,'xtick',[]);
        subplot(2,5,9);
        plot(dat{imarker}.time{1},shifted(rejected,:)');
        axis tight
        
        subplot(2,5,5);
        imagesc(shifted_clean_z);
        set(gca,'xtick',[]);
        subplot(2,5,10);
        plot(dat{imarker}.time{1},shifted_clean_z');
        axis tight
        
        set(fig,'PaperOrientation','landscape');
        set(fig,'PaperUnits','normalized');
        set(fig,'PaperPosition', [0 0 1 1]);
        print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'alignment_',cfg.name{imarker},'.pdf']),'-r600');
        print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,'alignment_',cfg.name{imarker},'.png']),'-r600');
        close all
        
        %% shift trialinfo in musestruct
        %             i = 1;
        %             hasmarker = zeros(length(MuseStruct_micro),1);
        %             for idir = 1:length(MuseStruct_micro)
        %                 if isfield(MuseStruct_micro{idir},'markers')
        %                     if isfield(MuseStruct_micro{idir}.markers,(cfg.muse.startend{imarker,1}))
        %                         if ~isempty(MuseStruct_micro{idir}.markers.(cfg.muse.startend{imarker,1}).events)
        %                             MuseStruct_micro{idir}.markers.(cfg.muse.startend{imarker,1}).events.sample = MuseStruct_micro{idir}.markers.(cfg.muse.startend{imarker,1}).events + nshift(i);
        %                             MuseStruct_micro{idir}.markers.(cfg.muse.startend{imarker,1}).events.offset = MuseStruct_micro{idir}.markers.(cfg.muse.startend{imarker,1}).offset + nshift(i);
        %                             MuseStruct_micro{idir}.markers.(cfg.muse.startend{imarker,1}).events.time   = MuseStruct_micro{idir}.markers.(cfg.muse.startend{imarker,1}).time   + dat_micro{1}.fsample / nshift(i);
        %                             i = i + 1;
        %                         end
        %                     end
        %                 end
        %             end
        
        
        
%         
%         i = 1;
%         clear q
%         hasmarker = zeros(length(MuseStruct_macro),1);
%         for idir = 1:length(MuseStruct_macro)
%             if isfield(MuseStruct_macro{idir},'markers')
%                 if isfield(MuseStruct_macro{idir}.markers,(cfg.muse.startend{imarker,1}))
%                     if ~isempty(MuseStruct_macro{idir}.markers.(cfg.muse.startend{imarker,1}).events)
%                         for ievent = 1 : size(MuseStruct_macro{idir}.markers.(cfg.muse.startend{imarker,1}).events,2)
% %                             MuseStruct_macro{idir}.markers.(cfg.muse.startend{imarker,1}).events(ievent).sample = MuseStruct_macro{idir}.markers.(cfg.muse.startend{imarker,1}).events(ievent).sample + nshift(i);
% %                             MuseStruct_macro{idir}.markers.(cfg.muse.startend{imarker,1}).events(ievent).offset = MuseStruct_macro{idir}.markers.(cfg.muse.startend{imarker,1}).events(ievent).offset + nshift(i);
% %                             MuseStruct_macro{idir}.markers.(cfg.muse.startend{imarker,1}).events(ievent).time   = MuseStruct_macro{idir}.markers.(cfg.muse.startend{imarker,1}).events(ievent).time   + dat_macro{1}.fsample / nshift(i);
%                                 
%                                    
%                                 i = i + 1;
%                                 q(i) = ievent;
%                         end
%                     end
%                 end
%             end
%         end
%         
        
        
        i = 1;
        clear q
        iidir = 1;
        hasmarker = zeros(length(MuseStruct_macro),1);
        for idir = 1:length(MuseStruct_macro)
            if isfield(MuseStruct_macro{idir},'markers')
                if isfield(MuseStruct_macro{idir}.markers,(cfg.muse.startend{imarker,1}))
                    if ~isempty(MuseStruct_macro{idir}.markers.(cfg.muse.startend{imarker,1}).events)
                        todelete = [];
                        for ievent = 1 : size(MuseStruct_macro{idir}.markers.(cfg.muse.startend{imarker,1}).events,2)
                            if any(dat{imarker}.trialinfo(:,2) == ievent & dat{imarker}.trialinfo(:,3) == idir    )
                                q(i) = ievent;
                                
                                MuseStruct_macro{idir}.markers.(cfg.muse.startend{imarker,1}).events(ievent).sample = MuseStruct_macro{idir}.markers.(cfg.muse.startend{imarker,1}).events(ievent).sample + nshift(i);
                                MuseStruct_macro{idir}.markers.(cfg.muse.startend{imarker,1}).events(ievent).offset = MuseStruct_macro{idir}.markers.(cfg.muse.startend{imarker,1}).events(ievent).offset + nshift(i);
                                MuseStruct_macro{idir}.markers.(cfg.muse.startend{imarker,1}).events(ievent).time   = MuseStruct_macro{idir}.markers.(cfg.muse.startend{imarker,1}).events(ievent).time   + dat_macro{1}.fsample / nshift(i);
                                i = i + 1;
                            else
                                todelete = [todelete, ievent];
                                
                            end
                        end
                         MuseStruct_macro{idir}.markers.(cfg.muse.startend{imarker,1}).events(todelete) = [];
                    end
                end
            end
        end
            
        
        save(fullfile(cfg.datasavedir,[cfg.prefix,'_MuseStruct_macro_aligned.mat']),'MuseStruct_macro');
        
        % load  and segmend data according to realiged timings
        [dat_micro, dat_macro] = readLFP(cfg, MuseStruct_micro, MuseStruct_macro, true, false);
        
        
        % baseline correct
        cfgtemp = [];
        cfgtemp.demean = 'yes';
        cfgtemp.baselinewindow = [-0.6 -0.3];
        dat{imarker} = ft_preprocessing(cfgtemp,dat_macro{imarker});
%         
%         % select data
%         cfgtemp = [];
%         cfgtemp.channel = 1;
%         cfgtemp.latency = [-0.1, 0.2];
%         dat{imarker}  = ft_selectdata(cfgtemp,dat{imarker});
        
        % put all in matrix
        clear input
        for itrial = 1 : size(dat{imarker}.trial,2)
            input(itrial,:) = dat{imarker}.trial{itrial}(1,:);
        end
        
        % z-normalize data
        z_input = nanznorm(input);
        z_input(isnan(z_input)) = 0;
        
        
        
        for N = 2 : 6
            
            %% Euclidian distance
            eucD = pdist(z_input,'seuclidean');
            %             eucD = pdist(shifted_clean_z,@naneucdist); % method to allow nans
%             
%             [Y,stress,disparities] = mdscale(eucD,3);
%             figure; 
%             scatter3(Y(:,1),Y(:,2),Y(:,3));
%             
            
            %% Hierarchical Clustering
            
            clustTreeEuc = linkage(eucD,'average');
            
            % plot tree
            fig = figure;
            fig.Renderer = 'Painters';
            
            subplot(1,3,1);
            [h,nodes] = dendrogram(clustTreeEuc,0);
            title(sprintf('Cophenetic corr: %0.2f',cophenet(clustTreeEuc,eucD)));
            set(gca,'xtick',[]);
            
            subplot(1,3,2);
            [h,nodes] = dendrogram(clustTreeEuc,10);
            title('Max leafs: 10');
            
            subplot(1,3,3);
            [h,nodes] = dendrogram(clustTreeEuc,N);
            title(sprintf('Max leafs: %d',N));
            xlabel('Leaf number');
            
            for i = 1 : N
                sumc(i) = sum(nodes==i);
            end
            set(gca,'Xticklabel',sumc);
            xlabel('#Datapoints');
            
            for i = 1 : N
                cluster{i}          = input(nodes==i,:);
                cluster_mean(i,:)   = mean(cluster{i});
            end
            ymin = min(min(cluster_mean));
            ymax = max(max(cluster_mean));
            
            set(fig,'PaperOrientation','landscape');
            set(fig,'PaperUnits','normalized');
            set(fig,'PaperPosition', [0 0 1 1]);
            print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'dendogram_',cfg.name{imarker},'_N',num2str(N),'.pdf']),'-r600');
            print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,'dendogram_',cfg.name{imarker},'_N',num2str(N),'.png']),'-r600');
            close all
            
            fig = figure;
            fig.Renderer = 'Painters';
            
            % plot grids
            for i = 1 : N
                
                [z_input_realigned,~] = align_data(z_input(nodes==i,:),5);
                
                subplot(3,N,i);
                imagesc(z_input(nodes==i,:));
                set(gca,'xtick',[]);
                
                subplot(3,N,i+N*1);
                imagesc(z_input_realigned);
                set(gca,'xtick',[]);
                
                subplot(3,N,i+N*2); hold;
                plot(dat{imarker}.time{1},nanmean(z_input(nodes==i,:)));
                %             xlim([dat{imarker}.time{1}(1),dat{imarker}.time{1}(end)]);
                %             ylim([ymin*1.1,ymax*1.1]);
                axis tight
                plot(dat{imarker}.time{1},nanmean(z_input_realigned));
                
            end
            %         subplot(3,N,N*2+1);
            %         [silh,h] = silhouette(shifted_clean_z,nodes,'euclidean');
            %
            set(fig,'PaperOrientation','landscape');
            set(fig,'PaperUnits','normalized');
            set(fig,'PaperPosition', [0 0 1 1]);
            print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'hierach_clusters_',cfg.name{imarker},'_N',num2str(N),'.pdf']),'-r600');
            print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,'hierach_clusters_',cfg.name{imarker},'_N',num2str(N),'.png']),'-r600');
            close all
            
            
            
            
            
            %% kmeans
            
            %         [indx,centroids]=kmeans(shifted_clean,N,'MaxIter',2000,'Replicates',100,'dist','sqeuclidean');
            [indx,centroids]=kmeans(z_input,N,'MaxIter',2000,'Replicates',100,'dist','sqeuclidean');
            
            clear cluster cluster_mean
            for i = 1 : N
                cluster{i}          = z_input(indx==i,:);
                cluster_mean(i,:)   = mean(cluster{i});
            end
            ymin = min(min(cluster_mean))*1.3;
            ymax = max(max(cluster_mean))*1.3;
            
            
            fig = figure; hold;
            fig.Renderer = 'Painters';
            
            for i = 1 : N
                
                [z_input_realigned,~] = align_data(z_input(indx==i,:),5);
                
                subplot(3,N,i);
                imagesc(z_input(indx==i,:));
                set(gca,'xtick',[]);
                
                subplot(3,N,i+N*1);
                imagesc(z_input_realigned);
                set(gca,'xtick',[]);
                
                subplot(3,N,i+N*2); hold;
                plot(dat{imarker}.time{1},nanmean(z_input(indx==i,:)));
                %             xlim([dat{imarker}.time{1}(1),dat{imarker}.time{1}(end)]);
                %             ylim([ymin*1.1,ymax*1.1]);
                axis tight
                plot(dat{imarker}.time{1},nanmean(z_input_realigned));
                
            end
            
            set(fig,'PaperOrientation','landscape');
            set(fig,'PaperUnits','normalized');
            set(fig,'PaperPosition', [0 0 1 1]);
            print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'raster_trials_kmeans',cfg.name{imarker},'_N',num2str(N),'.pdf']),'-r600');
            print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,'raster_trials_kmeans',cfg.name{imarker},'_N',num2str(N),'.png']),'-r600');
            close all
            
            % plot siluette
            fig = figure;
            [silh,h] = silhouette(z_input,indx,'sqeuclidean');
            set(fig,'PaperOrientation','landscape');
            set(fig,'PaperUnits','normalized');
            set(fig,'PaperPosition', [0 0 1 1]);
            print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'raster_trials_kmeans_silhouette',cfg.name{imarker},'_N',num2str(N),'.pdf']),'-r600');
            print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,'raster_trials_kmeans_silhouette',cfg.name{imarker},'_N',num2str(N),'.png']),'-r600');
            close all
            
            
            
            
            
            
            %% Dynamic Time Warping
            
            maxamp = 100;
            tic;
            
            if maxamp == 0
                dtwD = pdist(z_input,@(Xi,Xj) dtwdist(Xi,Xj,'squared'));
            else
                dtwD = pdist(z_input,@(Xi,Xj) dtwdist(Xi,Xj,maxamp,'squared'));
            end
            
            save(fullfile(cfg.datasavedir,[cfg.prefix,'dtwD_sel.mat']),'dtwD');
            elapsedTime = toc
            
            fig = figure;
            fig.Renderer = 'Painters';
            imagesc(squareform(dtwD));
            colorbar;
            title('Distance Matrix');
            set(fig,'PaperOrientation','landscape');
            set(fig,'PaperUnits','normalized');
            set(fig,'PaperPosition', [0 0 1 1]);
            print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'distance_matrix_DTW_',cfg.name{imarker},'_N',num2str(N),'.pdf']),'-r600');
            print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,'distance_matrix_DTW_',cfg.name{imarker},'_N',num2str(N),'.png']),'-r600');
            close all
            
            
            % plot tree
            
            clustTreeEuc = linkage(dtwD,'weighted');
            clear sumc
            
            fig = figure;
            fig.Renderer = 'Painters';
            
            subplot(1,3,1);
            [h,nodes] = dendrogram(clustTreeEuc,0);
            title(sprintf('Cophenetic corr: %0.2f',cophenet(clustTreeEuc,eucD)));
            set(gca,'xtick',[]);
            
            subplot(1,3,2);
            [h,nodes] = dendrogram(clustTreeEuc,10);
            title('Max leafs: 10');
            
            subplot(1,3,3);
            [h,nodes] = dendrogram(clustTreeEuc,N);
            title(sprintf('Max leafs: %d',N));
            xlabel('Leaf number');
            
            for i = 1 : N
                sumc(i) = sum(nodes==i);
            end
            set(gca,'Xticklabel',sumc);
            xlabel('#Datapoints');
            
            %         for i = 1 : N
            %             cluster{i}          = input(nodes==i,:);
            %             cluster_mean(i,:)   = mean(cluster{i});
            %         end
            %         ymin = min(min(cluster_mean));
            %         ymax = max(max(cluster_mean));
            %
            
            set(fig,'PaperOrientation','landscape');
            set(fig,'PaperUnits','normalized');
            set(fig,'PaperPosition', [0 0 1 1]);
            print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'dendogram_DTW_',cfg.name{imarker},'.pdf']),'-r600');
            print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,'dendogram_DTW_',cfg.name{imarker},'.png']),'-r600');
            close all
            
            % plot grids
            fig = figure;
            fig.Renderer = 'Painters';
            
            for i = 1 : N
                
                [z_input_realigned,~] = align_data(z_input(indx==i,:),5);
                
                subplot(3,N,i);
                imagesc(z_input(indx==i,:));
                set(gca,'xtick',[]);
                
                subplot(3,N,i+N*1);
                imagesc(z_input_realigned);
                set(gca,'xtick',[]);
                
                subplot(3,N,i+N*2); hold;
                plot(dat{imarker}.time{1},nanmean(z_input(indx==i,:)));
                %             xlim([dat{imarker}.time{1}(1),dat{imarker}.time{1}(end)]);
                %             ylim([ymin*1.1,ymax*1.1]);
                axis tight
                plot(dat{imarker}.time{1},nanmean(z_input_realigned));
                
            end
            %         subplot(3,N,N*2+1);
            %         [silh,h] = silhouette(shifted_clean_z,nodes,'euclidean');
            
            set(fig,'PaperOrientation','landscape');
            set(fig,'PaperUnits','normalized');
            set(fig,'PaperPosition', [0 0 1 1]);
            print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'hierach_clusters_DTW_',cfg.name{imarker},'_N',num2str(N),'.pdf']),'-r600');
            print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,'hierach_clusters_DTW_',cfg.name{imarker},'_N',num2str(N),'.png']),'-r600');
            close all
            
        end
        
    end
    
end