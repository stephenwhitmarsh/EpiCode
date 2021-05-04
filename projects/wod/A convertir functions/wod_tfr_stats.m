


% %% Make comparison between electrodes
% 
% %arrange data for ANOVA
% for idata= 1:size(analysis_names,2)
%     for ichan= 1:size(chan_list,1)
%         chan_name= chan_list{ichan};
%         
%         for iband =[ "HF" "MF" "MLF" "LF"]
%             
%             %for timing
%             stat.(iband).(analysis_names{idata}).time(:,ichan)= time.(iband).(analysis_names{idata}).(chan_name);
%             
%             %for values
%             stat.(iband).(analysis_names{idata}).value(:,ichan)= value.(iband).(analysis_names{idata}).(chan_name);
%             
%         end %iband
%     end %ichan
%     
%     for iband =[ "HF" "MF" "MLF" "LF"]
%         
%         %Express timings with earliest peak at 0
%         for irat= 1:size(stat.(iband).(analysis_names{idata}).time,1)
%             stat.(iband).(analysis_names{idata}).time_norm(irat,:)=  stat.(iband).(analysis_names{idata}).time(irat,:)-min(stat.(iband).(analysis_names{idata}).time(irat,:));
%         end %irat
%         
%         %for timings
%         [p, tbl,stats]= kruskalwallis(stat.(iband).(analysis_names{idata}).time_norm,[],'off');
%         %for values
%         [p1, tbl1, stats1]= kruskalwallis(stat.(iband).(analysis_names{idata}).value,[],'off');
%         
%         %make multiple comparison
%         %for timings
%         [c, m]= multcompare(stats);
%         
%         %for values
%         [c1, m1]= multcompare(stats1);
%         
%         %store data in structure
%         Stats_depth.(analysis_names{idata}).time.(iband).p= p;
%         Stats_depth.(analysis_names{idata}).time.(iband).comp(:,:)= c;
%         Stats_depth.(analysis_names{idata}).time.(iband).mean(:,:)= m;
%         
%         Stats_depth.(analysis_names{idata}).value.(iband).p= p1;
%         Stats_depth.(analysis_names{idata}).value.(iband).comp(:,:)= c1;
%         Stats_depth.(analysis_names{idata}).value.(iband).mean(:,:)= m1;
%         
%         clear p1 p stats1 stats c c1
%         
%     end %iband
% end %idata
% 
% save(fullfile(freqstatpath,'Stats_depth.mat'),'Stats_depth');
% save(fullfile(freqstatpath,'all_data_norm.mat'),'stat');
% clear Stats_depth time value

%% Make Spearman corr of timings with depth

% %load depth into structure (glisser table dans variables)
% Detectionpath=fullfile(config{4}.datasavedir,'Detection');
% temp= load(fullfile(Detectionpath,'Depth_electrode.mat'));
% Depth= temp.Electrode_depth;
% clear temp
% 
% a=0
% for irat=4:size(config,2)
%     
%     a=a+1;
%     if irat== size(config,2)
%         continue
%     end
%     A=Depth.(sprintf('Rat_%i',irat));
%     B=Depth.(sprintf('Rat_%i',irat+1));
%     C= horzcat(A,B)
%     
% end %irat
% 
% for idata= 1: size(analysis_names,2)
%     for iband= ["HF" "MF" "MLF" "LF"]
%         for ival= ["time" "time_norm"]
%             %delete last electrode to have same number
%             stat.(iband).(analysis_names{idata}).(ival)(:,end)=[];
%             
%             %Make Spearman corr
%             [rho,pval]= corr(transpose(stat.depth),transpose(stat.(iband).(analysis_names{idata}).(ival)),'Type','Spearman','Rows','complete');
%             Corr.Rho.(analysis_names{idata}).(ival).(iband).all=rho;
%             Corr.Pval.(analysis_names{idata}).(ival).(iband).all=pval;
%             
%             clear rho pval
%             
%             if ~isfolder(fullfile(statsavedir,'corr_depth'))
%                 mkdir(fullfile(statsavedir,'corr_depth'))
%             end
%             
%             %save Rho and P-values for each protocol in excel files
%             writematrix(Corr.Rho.(analysis_names{idata}).(ival).(iband).all(1,:),fullfile(statsavedir,'corr_depth',sprintf('%s_%s_%s_Rho_depth.csv',analysis_names{idata},(ival),(iband))));
%             writematrix(Corr.Pval.(analysis_names{idata}).(ival).(iband).all(1,:),fullfile(statsavedir,'corr_depth',sprintf('%s_%s_%s_Pval_depth.csv',analysis_names{idata},(ival),(iband))));
%             
%             %Averaging Corr values and saving structure
%             Corr.Rho.(analysis_names{idata}).(ival).(iband).mean(1,1)=nanmean(Corr.Rho.(analysis_names{idata}).(ival).(iband).all(1,:));
%             Corr.Pval.(analysis_names{idata}).(ival).(iband).mean(1,1)=nanmean(Corr.Pval.(analysis_names{idata}).(ival).(iband).all(1,:));
%             save(fullfile(statsavedir,'corr_depth','all_Corr.mat'),'Corr');
%             
%         end %ival
%     end %iband
% end %idata
% 
% %% Compare WoD delay with peak times delay
% 
% %load WoD delay into structure (glisser table dans variables)
% A= table2array(woddelay);
% A(1,:)=[];
% A=transpose(A);
% for ichan= 1: size(stat.HF.timefreq_wod.time,2)
%     stat.wod(:,ichan)=A(:,ichan);
% end
% clear A
% % Attribute to each peak time a WoD delay
% for idata= 1:size(analysis_names,2)
%     for iband= ["HF" "MF" "MLF" "LF"]
%         for ival= ["time" "time_norm"]
%             
%             A=stat.wod(:);
%             B= stat.(iband).(analysis_names{idata}).(ival)(:);
%             Comp_wod.(iband).(analysis_names{idata}).(ival)=[A B];
%             clear A B
%             
%         end %ival
%     end %iband
% end %idata
% 
% %Save structure to plot
% 
% comp_woddir=fullfile(statsavedir,'comp_wod');
% 
% if ~isfolder(comp_woddir)
%     mkdir(comp_woddir);
% end
% 
% save(fullfile(comp_woddir,'comp_wod_toplot.mat'),'Comp_wod');
% 
% %Make Spearman Corr with wod delay
% 
% stat.wod(:,end)=[];
% 
% 
% for idata= 1:size(analysis_names,2)
%     for iband= ["HF" "MF" "MLF" "LF"]
%         for ival= ["time" "time_norm"]
%             %stat.(iband).(analysis_names{idata}).(ival)(:,end)=[];
%             
%             %Make Spearman corr
%             [rho,pval]= corr(transpose(stat.wod),transpose(stat.(iband).(analysis_names{idata}).(ival)),'Type','Spearman','Rows','complete');
%             Corr_wod.Rho.(analysis_names{idata}).(ival).(iband).all=rho;
%             Corr_wod.Pval.(analysis_names{idata}).(ival).(iband).all=pval;
%             
%             Corr_wod.Rho.(analysis_names{idata}).(ival).(iband).mean=nanmean(Corr_wod.Rho.(analysis_names{idata}).(ival).(iband).all);
%             Corr_wod.Pval.(analysis_names{idata}).(ival).(iband).mean=nanmean(Corr_wod.Pval.(analysis_names{idata}).(ival).(iband).all);
%             
%             
%         end %ival
%     end %iband
% end %idata
% 
% %save data
% save(fullfile(comp_woddir,'corr_wod.mat'),'Corr_wod');

%% Compare delays between frequency bands
bands={"HF", "MF", "MLF", "LF"};

%between adjacent bands
for idata= 1:size(analysis_names,2)
    for iband= 1:size(bands,2)
        for ichan= 1:size(stat.HF.timefreq_wod.time,2)
            
            if iband >= size(bands,2)
                
                continue
                
            end
            
            A=stat.(bands{iband}).(analysis_names{idata}).time(:,ichan);
            B=stat.(bands{iband+1}).(analysis_names{idata}).time(:,ichan);
            
            [p,h,stats]= signrank(A,B);
            
            Test_betbands.(bands{iband}).(analysis_names{idata}).p(ichan,1)= p;
            Test_betbands.(bands{iband}).(analysis_names{idata}).h(ichan,1)=h;
            Test_betbands.(bands{iband}).(analysis_names{idata}).stats=stats;
            
            clear A B p h stats
            
        end %ichan
        if iband >= size(bands,2)
                
                continue
                
            end
        
        [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(Test_betbands.(bands{iband}).(analysis_names{idata}).p,0.05,'pdep','no');
        Test_betbands.(bands{iband}).(analysis_names{idata}).correc_p= adj_p;
        
    end %iband
end %idata

%between ends of the spectrum
bands={"HF","LF"};

for idata= 1:size(analysis_names,2)
    for iband= 1:size(bands,2)
        for ichan= 1:size(stat.HF.timefreq_wod.time,2)
            
            if iband >= size(bands,2)
                continue
            end
            
            A=stat.(bands{iband}).(analysis_names{idata}).time(:,ichan);
            B=stat.(bands{iband+1}).(analysis_names{idata}).time(:,ichan);
            
            [p,h,stats]= signrank(A,B);
            
            Test_betbands.extreme.(bands{iband}).(analysis_names{idata}).p(ichan,1)= p;
            Test_betbands.extreme.(bands{iband}).(analysis_names{idata}).h(ichan,1)=h;
            Test_betbands.extreme.(bands{iband}).(analysis_names{idata}).stats=stats;
            
            clear A B p h stats
            
        end %ichan
        
        %correction of p-values
        [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(Test_betbands.extreme.(bands{iband}).(analysis_names{idata}).p,0.05,'pdep','no');
        Test_betbands.extreme.(bands{iband}).(analysis_names{idata}).correc_p= adj_p;
    end %iband
end %idata
clear bands

save(fullfile(freqstatpath,'test_between_freqband.mat'),'Test_betbands');