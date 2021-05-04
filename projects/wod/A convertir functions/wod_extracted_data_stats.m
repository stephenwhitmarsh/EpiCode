function statistics_propag_origin= wod_extracted_data_stats(cfg,force)


detect_stats_dir= fullfile(cfg{4}.statsavedir,'Waves_detections');
fname_out=fullfile(detect_stats_dir,'calculateddata_stats.mat');

if exist(fname_out, 'file') && force == false
    load(fname_out, 'calculateddata_stats');
    return
end

%% Load and pool calculated data

detectiondatapath= fullfile(cfg{4}.datasavedir,'Detection');
load(fullfile(detectiondatapath,'calculated_data.mat'),'Data_waves');

for irat= 1:size(cfg,2)
    
    if isempty(cfg{irat})
        continue
    end
    
    for iwave= ["WOD" "WOR"]
        for iana= ["peak_time" "min_slope" "start_time"]
            
            first_trial= Data_waves{irat}.(iwave).(iana)
            
            
        end %iana
    end %iwave
end %irat





%% Make statiscal analysis between trials

% origin time and depth
for iwave=["WOD" "WOR"]
    for iana= ["origin_time" "origin_depth"]
        count=0;
        for itime= ["peak_time" "min_slope" "start_time"]
            count= count+1;
            %separate trials as different arrays
            first_trial=Data_waves.(iwave).(iana).(itime)(:,1);
            second_trial=Data_waves.(iwave).(iana).(itime)(:,2);
            %t-test between trials
            p=signrank(first_trial,second_trial);
            p_val_trials.(iwave).(iana)(count,1)=p;
        end %itime
        %p values corrections
        [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p_val_trials.(iwave).(iana));
        adj_pval_trials.(iwave).(iana)=adj_p;
    end %iana
end %iwave
clear first_trial second_trial

%same operation for average propagation speed
for iwave=["WOD" "WOR"]
    for itime= ["peak_time" "min_slope" "start_time"]
        count=0;
        
        for isens=["up" "down"]
            
            count=count+1;
            %separate trials as different arrays
            first_trial=Data_waves.(iwave).speed.(itime).(isens)(:,1);
            second_trial=Data_waves.(iwave).speed.(itime).(isens)(:,2);
            %t-test between trials
            p=ranksum(first_trial,second_trial);
            p_val_trials.(iwave).speed.(itime).(isens)(count,1)=p;
            
        end %isens
    end %itime
end %iwave


%% Make statistical analysis between waves

for iana= ["origin_time" "origin_depth"]
    count=0;
    for itime= ["peak" "min_slope" "start"]
        count= count+1;
        wod_first_trial=Data_waves.WOD.(iana).(itime)(:,1);
        wod_second_trial=Data_waves.WOD.(iana).(itime)(:,2);
        wod_alltrials=vertcat(wod_first_trial,wod_second_trial);
        
        wor_first_trial=Data_waves.WOR.(iana).(itime)(:,1);
        wor_second_trial=Data_waves.WOR.(iana).(itime)(:,2);
        wor_alltrials=vertcat(wor_first_trial,wor_second_trial);
        
        p=signrank(wod_alltrials,wor_alltrials);
        p_val_waves.(iana)(count,1)=p;
        
    end %itime
    [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p_val_waves.(iana));
    adj_pval_waves.(iana)=adj_p;
end %iana

detect_stats_dir= fullfile(cfg{4}.statsavedir,'Waves_detections');

if ~isfolder(detect_stats_dir)
    mkdir(detect_stats_dir);
end

%store p-values in a single structure
Stats.pval_waves=p_val_waves;
Stats.adj_pval_waves=adj_pval_waves;

Stats.pval_trials=p_val_trials;
Stats.adj_pval_trials=adj_pval_trials;


save(fullfile(detect_stats_dir,'tests_mean_calculated_data.mat'),'Stats');

end %function