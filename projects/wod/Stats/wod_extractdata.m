
if ispc
    addpath (genpath('\\lexport\iss01.charpier\analyses\wod\Antoine\EpiCode\shared'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\wod\Antoine\EpiCode\external'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\wod\Antoine\EpiCode\templates'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\wod\Antoine\EpiCode\projects\wod'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\wod\Antoine\EpiCode\projects\dtx'))
    addpath \\lexport\iss01.charpier\analyses\wod\fieldtrip-20200607
    addpath \\lexport\iss01.charpier\analyses\wod\fdr_bh
elseif isunix
    addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/Antoine/EpiCode/shared'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/Antoine/EpiCode/external'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/Antoine/EpiCode/templates'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/Antoine/EpiCode/projects/wod'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/Antoine/EpiCode/projects/dtx'))
    addpath /network/lustre/iss01/charpier/analyses/wod/fieldtrip-20200607
end

ft_defaults

config = wod_setparams;

freqstatpath= fullfile(config{4}.statsavedir,'freq_data');
analysis_names={'timefreq_wod','timefreq_wod_timenorm','timefreq_wod_blcorrected','timefreq_wod_timenorm_blcorrected'};

%% Prepare structure to analyze

for idata= 1: size(analysis_names,2)
    for irat= 4:size(config,2)
        cfgcommon = config{4}; %same common config for all rats
        
        fprintf('Load timefreq data for rat %d/%d\n', irat,size(config,2));
        data_temp = load(fullfile(config{irat}.datasavedir,[config{irat}.prefix,analysis_names{idata},'.mat']));
        
        data_rat = [];
        count_trials                    = 0;
        %get data from one rat in a structure
        for itrial = 1:size(data_temp.(analysis_names{idata}),2)
            count_trials = count_trials +1;
            wod_rat(count_trials) = count_trials; %store rat id for each wod (some rats have several wod)
            chan_list                                   = fieldnames(data_temp.(analysis_names{idata}){itrial});
            
            data_rat= data_temp.(analysis_names{idata}){itrial};
            
            for iband= ["HF" "MF" "MLF" "LF"]
                for ichan= 1: numel(chan_list)
                    chan_name = chan_list{ichan};
                    %replace empty channels by nans
                    hasdata = find(~structfun(@isempty,data_rat));
                    if isempty(hasdata)
                        continue
                    end
                    
                    if isempty(data_rat.(chan_name))
                        data_rat.(chan_name)                       = data_rat.(chan_list{hasdata(1)});
                        data_rat.(chan_name).powspctrm             = ones(size(data_rat.(chan_name).powspctrm));
                    end
                    
                    
                    %% Select frequency bands
                    
                    
                    
                    % Select frequencies and define time window to search
                    % peak
                    
                    t1=data_rat.(chan_name).time(1)+10;
                    t2=data_rat.(chan_name).time(end)-30;
                    tsearch= [t1 t2];
                    
                    cfgtemp= [];
                    cfgtemp.latency= tsearch;
                    cfgtemp.frequency= config{irat}.timefreq.(iband);
                    %cfgtemp.avgoverfreq = 'yes';
                    pow.(iband).(chan_name)= ft_selectdata(cfgtemp, data_rat.(chan_name));
                    %average over frequency
                    pow.(iband).(chan_name).powspctrm = permute(pow.(iband).(chan_name).powspctrm, [2 3 1]); %put first dimension as last so it is deleted
                    pow.(iband).(chan_name).powspctrm = nanmean(pow.(iband).(chan_name).powspctrm,1);
                    pow.(iband).(chan_name).time= pow.(iband).(chan_name).time;
                    
                    while iband= "LF"
                    Ratio.(chan_name)= pow.(iband).(chan_name);
                    Ratio.(chan_name).powspctrm= pow.(iband).(chan_name).powspctrm ./ pow.HF.(chan_name).powspctrm;
                    end
                    
                    
                    %% find peak value and peak time
                    
                    
                    [v_peak.(iband).(chan_name), t_peak.(iband).(chan_name)] = findpeaks(pow.(iband).(chan_name).powspctrm,pow.(iband).(chan_name).time,'NPeaks',1,'SortStr','descend','WidthReference','Halfheight');
                    if isempty([v_peak.(iband).(chan_name), t_peak.(iband).(chan_name)])
                        v_peak.(iband).(chan_name)= NaN;
                        t_peak.(iband).(chan_name)= NaN;
                    end
                    
                    %plot detection of peaks
                    %                     fig= figure;
                    %                     plot(pow.(iband).(chan_name).time,pow.(iband).(chan_name).powspctrm)
                    %                     hold on
                    %                     scatter(t_peak
                    
                    
                    %store values in structure
                    Val.(analysis_names{idata}).(chan_name).(iband).peak_value(irat,itrial)= v_peak.(iband).(chan_name);
                    Val.(analysis_names{idata}).(chan_name).(iband).peak_time(irat,itrial)= t_peak.(iband).(chan_name);
                    
                    %% Rearrange data for statistic analysis
                    
                    %for peak time and values
                    time.(iband).(analysis_names{idata}).(chan_name)= nonzeros(Val.(analysis_names{idata}).(chan_name).(iband).peak_time);
                    value.(iband).(analysis_names{idata}).(chan_name)= nonzeros(Val.(analysis_names{idata}).(chan_name).(iband).peak_value);
                    
                    
                    
                    
                end %ichan
                
                %un plot avec toutes les électrodes de la bande + toutes
                %les détections (scatter)
                %plot detection of peaks
                
                %FIXME Save Figures in a detection folder
                fig= figure;
                for ichan= 1:numel(chan_list)
                    
                    chan_name= chan_list{ichan};
                    
                    
                    plot(pow.(iband).(chan_name).time,pow.(iband).(chan_name).powspctrm);
                    hold on
                    scatter(t_peak.(iband).(chan_name),v_peak.(iband).(chan_name),'rx');
                    xlim([0 60])
                end
                
                %save figures
                fname= fullfile(config{4}.imagesavedir_data{2},'detection',sprintf('%sWOD%i_%s_%s',config{irat}.prefix,itrial,iband,analysis_names{idata}));
                dtx_savefigure(fig,fname,'png','pdf','close');
            end %iband
        end %itrial
    end %irat
end %idata

%% Save data

if ~isfolder(freqstatpath)
    mkdir(freqstatpath);
end

save(fullfile(freqstatpath,'peak_time.mat'),'time');
save(fullfile(freqstatpath,'peak_value.mat'),'value');
save(fullfile(freqstatpath,'all_data.mat'),'Val');




%% Make comparison between electrodes

%arrange data for ANOVA
for idata= 1:size(analysis_names,2)
    for ichan= 1:size(chan_list,1)
        chan_name= chan_list{ichan};
        
        
        for iband =[ "HF" "MF" "MLF" "LF"]
            
            %for timing
            stat.(iband).(analysis_names{idata}).time(:,ichan)= time.(iband).(analysis_names{idata}).(chan_name);
            
            %for values
            stat.(iband).(analysis_names{idata}).value(:,ichan)= value.(iband).(analysis_names{idata}).(chan_name);
            
        end %iband
    end %ichan
    
    for iband =[ "HF" "MF" "MLF" "LF"]
        
        %Express timings with earliest peak at 0
        for irat= 1:size(stat.(iband).(analysis_names{idata}).time,1)
            stat.(iband).(analysis_names{idata}).time_norm(irat,:)=  stat.(iband).(analysis_names{idata}).time(irat,:)-min(stat.(iband).(analysis_names{idata}).time(irat,:));
        end %irat
        
        
        
        %for timings
        [p, tbl,stats]= kruskalwallis(stat.(iband).(analysis_names{idata}).time_norm,[],'off');
        %for values
        [p1, tbl1, stats1]= kruskalwallis(stat.(iband).(analysis_names{idata}).value,[],'off');
        
        %make multiple comparison
        %for timings
        [c, m]= multcompare(stats);
        
        %for values
        [c1, m1]= multcompare(stats1);
        
        
        %store data in structure
        Stats_depth.(analysis_names{idata}).time.(iband).p= p;
        Stats_depth.(analysis_names{idata}).time.(iband).comp(:,:)= c;
        Stats_depth.(analysis_names{idata}).time.(iband).mean(:,:)= m;
        
        Stats_depth.(analysis_names{idata}).value.(iband).p= p1;
        Stats_depth.(analysis_names{idata}).value.(iband).comp(:,:)= c1;
        Stats_depth.(analysis_names{idata}).value.(iband).mean(:,:)= m1;
        
        clear p1 p stats1 stats c c1
        
    end %iband
end %idata

save(fullfile(freqstatpath,'Stats_depth.mat'),'Stats_depth');
save(fullfile(freqstatpath,'all_data_norm.mat'),'stat');
clear Stats_depth time value

%% Make Spearman corr of timings with depth

%load depth into structure (glisser table dans variables)
Detectionpath=fullfile(config{4}.datasavedir,'Detection');
temp= load(fullfile(Detectionpath,'Depth_electrode.mat'));
Depth= temp.Electrode_depth;
clear temp

a=0
for irat=4:size(config,2)
    
    a=a+1;
    if irat== size(config,2)
        continue
    end
    A=Depth.(sprintf('Rat_%i',irat));
    B=Depth.(sprintf('Rat_%i',irat+1));
    C= horzcat(A,B)
    
end %irat

for idata= 1: size(analysis_names,2)
    for iband= ["HF" "MF" "MLF" "LF"]
        for ival= ["time" "time_norm"]
            %delete last electrode to have same number
            stat.(iband).(analysis_names{idata}).(ival)(:,end)=[];
            
            %Make Spearman corr
            [rho,pval]= corr(transpose(stat.depth),transpose(stat.(iband).(analysis_names{idata}).(ival)),'Type','Spearman','Rows','complete');
            Corr.Rho.(analysis_names{idata}).(ival).(iband).all=rho;
            Corr.Pval.(analysis_names{idata}).(ival).(iband).all=pval;
            
            clear rho pval
            
            if ~isfolder(fullfile(statsavedir,'corr_depth'))
                mkdir(fullfile(statsavedir,'corr_depth'))
            end
            
            %save Rho and P-values for each protocol in excel files
            writematrix(Corr.Rho.(analysis_names{idata}).(ival).(iband).all(1,:),fullfile(statsavedir,'corr_depth',sprintf('%s_%s_%s_Rho_depth.csv',analysis_names{idata},(ival),(iband))));
            writematrix(Corr.Pval.(analysis_names{idata}).(ival).(iband).all(1,:),fullfile(statsavedir,'corr_depth',sprintf('%s_%s_%s_Pval_depth.csv',analysis_names{idata},(ival),(iband))));
            
            %Averaging Corr values and saving structure
            Corr.Rho.(analysis_names{idata}).(ival).(iband).mean(1,1)=nanmean(Corr.Rho.(analysis_names{idata}).(ival).(iband).all(1,:));
            Corr.Pval.(analysis_names{idata}).(ival).(iband).mean(1,1)=nanmean(Corr.Pval.(analysis_names{idata}).(ival).(iband).all(1,:));
            save(fullfile(statsavedir,'corr_depth','all_Corr.mat'),'Corr');
            
        end %ival
    end %iband
end %idata

%% Compare WoD delay with peak times delay

%load WoD delay into structure (glisser table dans variables)
A= table2array(woddelay);
A(1,:)=[];
A=transpose(A);
for ichan= 1: size(stat.HF.timefreq_wod.time,2)
    stat.wod(:,ichan)=A(:,ichan);
end
clear A
% Attribute to each peak time a WoD delay
for idata= 1:size(analysis_names,2)
    for iband= ["HF" "MF" "MLF" "LF"]
        for ival= ["time" "time_norm"]
            
            A=stat.wod(:);
            B= stat.(iband).(analysis_names{idata}).(ival)(:);
            Comp_wod.(iband).(analysis_names{idata}).(ival)=[A B];
            clear A B
            
        end %ival
    end %iband
end %idata

%Save structure to plot

comp_woddir=fullfile(statsavedir,'comp_wod');

if ~isfolder(comp_woddir)
    mkdir(comp_woddir);
end

save(fullfile(comp_woddir,'comp_wod_toplot.mat'),'Comp_wod');

%Make Spearman Corr with wod delay

stat.wod(:,end)=[];


for idata= 1:size(analysis_names,2)
    for iband= ["HF" "MF" "MLF" "LF"]
        for ival= ["time" "time_norm"]
            %stat.(iband).(analysis_names{idata}).(ival)(:,end)=[];
            
            %Make Spearman corr
            [rho,pval]= corr(transpose(stat.wod),transpose(stat.(iband).(analysis_names{idata}).(ival)),'Type','Spearman','Rows','complete');
            Corr_wod.Rho.(analysis_names{idata}).(ival).(iband).all=rho;
            Corr_wod.Pval.(analysis_names{idata}).(ival).(iband).all=pval;
            
            Corr_wod.Rho.(analysis_names{idata}).(ival).(iband).mean=nanmean(Corr_wod.Rho.(analysis_names{idata}).(ival).(iband).all);
            Corr_wod.Pval.(analysis_names{idata}).(ival).(iband).mean=nanmean(Corr_wod.Pval.(analysis_names{idata}).(ival).(iband).all);
            
            
        end %ival
    end %iband
end %idata

%save data
save(fullfile(comp_woddir,'corr_wod.mat'),'Corr_wod');

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
