function wod_propag_analysis(rat_list,configscript)

try %en local
    scriptpath = matlab.desktop.editor.getActiveFilename;
catch %cluster
    scriptpath = mfilename('fullpath');
end

epicodepath = [fileparts(fileparts(fileparts(scriptpath))), filesep];

addpath (genpath([epicodepath,'development']))
addpath (genpath([epicodepath,'shared']))
addpath (genpath([epicodepath,'external']))
addpath (genpath([epicodepath,'templates']))
addpath (genpath([epicodepath,'projects', filesep, 'wod']))
addpath (genpath([epicodepath,'projects', filesep, 'dtx']))

if ispc
    addpath \\lexport\iss01.charpier\analyses\wod\fieldtrip-20200607
elseif isunix
    addpath /network/lustre/iss01/charpier/analyses/wod/fieldtrip-20200607
end

ft_defaults

config = eval(configscript);
detectiondatapath= fullfile(config{4}.datasavedir,'Detection');


%% loading data
temp= load(fullfile(detectiondatapath,'WoD_data.mat'));
WoD_data=temp.WOD_data;
clear temp
temp= load(fullfile(detectiondatapath,'WoR_data.mat'));
WoR_data=temp.WOR_data;
clear temp
temp= load(fullfile(detectiondatapath,'Depth_electrode.mat'));
Depth=temp.Electrode_depth;
clear temp

for irat= 4:size(config,2)
    irat_name=sprintf('Rat_%i',irat);
    for itime= ["peak" "min_slope" "start"]
        for itrial= 1:size(WoD_data.timings.(irat_name).(itime),2)
            %% WOD Find origin and store origin depth and timings
            A=WoD_data.timings.(irat_name).(itime)(:,itrial);
            B= Depth.(irat_name)(:,itrial);
            %find minimum timing
            origin_time=min(A);
            %store minimum timing in structure
            Data_waves.WOD.origin_time.(itime)(irat,itrial)=origin_time;
            
            %find index of origin and get origin depth
            idx_origin=find(A==origin_time);
            origin_depth=Depth.(irat_name)(idx_origin,itrial);
            %store origin depth
            Data_waves.WOD.origin_depth.(itime)(irat,itrial)=origin_depth;
            
            %% WOD Calculate propagation speed
            
            %Calculate propagation speed from origin in both ways
            for ichan= 1:size(A,1)
                
                if ichan== idx_origin
                    continue
                end
                
                if ichan<idx_origin
                    Speed_up(ichan,1)= (B(idx_origin,1)-B(ichan,1))/(A(ichan,1)-A(idx_origin,1));
                end
                
                
                if ichan > idx_origin
                    Speed_down(ichan,1)= (B(ichan,1)-B(idx_origin,1))/(A(ichan,1)-A(idx_origin,1));
                end
            end %ichan
            
            %replace zeros by NaN
            exist Speed_up
            bool_up= ans;
            exist Speed_down
            bool_down= ans;
            if bool_up==1
                Nul_values_up= find(Speed_up(:,1)==0);
                Speed_up(Nul_values_up,1)=NaN;
            end
            
            if bool_down==1
                Nul_values_down= find(Speed_down(:,1)==0);
                Speed_down(Nul_values_down,1)=NaN;
            end
            %Average and store values
            if bool_up==1
                Data_waves.WOD.speed.(itime).up(irat,itrial)=nanmean(Speed_up);
            end
            if bool_down==1
                Data_waves.WOD.speed.(itime).down(irat,itrial)=nanmean(Speed_down);
            end
            
            %% WOD Calculate Instantaneous speed
            
            %upwards
            for ichan= 1:idx_origin
                
                if ichan== idx_origin
                    continue
                end
                Speed_instan(ichan,1)= (B(ichan+1,1)-B(ichan,1))/(A(ichan,1)-A(ichan+1,1));
                
            end %ichan
            
            %downwards
            
            for ichan= idx_origin:size(A,1)
                
                if ichan==size(A,1)
                continue
                end
                
                Speed_instan(ichan,1)= (B(ichan+1,1)-B(ichan,1))/(A(ichan+1,1)-A(ichan,1));
               
            end %ichan
            
            Data_waves.WOD.speed_instan.(sprintf('Rat_%i',irat)).(itime)(:,itrial)=Speed_instan;
            
            
            
            clear A origin_time origin_depth bool_up bool_down Speed_up Speed_down Speed_instan
            %% WOR Find origin and store origin depth and timings
            
            %same operation for WoR
            A=WoR_data.timings.(irat_name).(itime)(:,itrial);
            
            %find minimum timing
            origin_time=min(A);
            %store minimum timing in structure
            Data_waves.WOR.origin_time.(itime)(irat,itrial)=origin_time;
            
            %find index of origin and get origin depth
            idx_origin=find(A==origin_time);
            origin_depth=Depth.(irat_name)(idx_origin,itrial);
            %store origin depth
            Data_waves.WOR.origin_depth.(itime)(irat,itrial)=origin_depth;
            
            %% WOR Calculate propagation speed
            
            %Calculate propagation speed from origin in both ways
            for ichan= 1:size(A,1)
                
                if ichan== idx_origin
                    continue
                end
                
                if ichan<idx_origin
                    Speed_up(ichan,1)= (B(idx_origin,1)-B(ichan,1))/(A(ichan,1)-A(idx_origin,1));
                end
                
                if ichan > idx_origin
                    Speed_down(ichan,1)= (B(ichan,1)-B(idx_origin,1))/(A(ichan,1)-A(idx_origin,1));
                end
            end %ichan
            
            %replace zeros by NaN
            exist Speed_up;
            bool_up= ans;
            exist Speed_down;
            bool_down= ans;
            if bool_up==1
                Nul_values_up= find(Speed_up==0);
                Speed_up(Nul_values_up,1)=NaN;
            end
            
            if bool_down==1
                Nul_values_down= find(Speed_down==0);
                Speed_down(Nul_values_down,1)=NaN;
            end
            %Average and store values
            if bool_up==1
                Data_waves.WOR.speed.(itime).up(irat,itrial)=nanmean(Speed_up);
            end
            if bool_down==1
                Data_waves.WOR.speed.(itime).down(irat,itrial)=nanmean(Speed_down);
            end
            
            %% WOR Calculate Instantaneous speed
            
            %upwards
            for ichan= 1:idx_origin
                
                if ichan== idx_origin
                    continue
                end
                Speed_instan(ichan,1)= (B(ichan+1,1)-B(ichan,1))/(A(ichan,1)-A(ichan+1,1));
                
            end %ichan
            
            %downwards
            
            for ichan= idx_origin:size(A,1)
                
                if ichan==size(A,1)
                continue
                end
                
                Speed_instan(ichan,1)= (B(ichan+1,1)-B(ichan,1))/(A(ichan+1,1)-A(ichan,1));
               
            end %ichan
            
            Data_waves.WOR.speed_instan.(sprintf('Rat_%i',irat)).(itime)(:,itrial)=Speed_instan;
            
            
            clear A origin_time origin_depth idx_origin bool_up bool_down Speed_up Speed_down Speed_instan
       
        end %itrial
    end %itime
end %irat

%% Replace all 0 values by NaNs
for iwave=["WOD" "WOR"]
    for iana= ["origin_time" "origin_depth"]
        for itime= ["peak" "min_slope" "start"]
            idx=find(Data_waves.(iwave).(iana).(itime)(:,:)==0);
            Data_waves.(iwave).(iana).(itime)(idx)=NaN;
        end %itime
    end %iana
end %iwave

save(fullfile(detectiondatapath,'calculated_data.mat'),'Data_waves');

%% Make statiscal analysis between trials


% origin time and depth
for iwave=["WOD" "WOR"]
    for iana= ["origin_time" "origin_depth"]
        count=0;
        for itime= ["peak" "min_slope" "start"]
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
    for itime= ["peak" "min_slope" "start"]
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

detect_stats_dir= fullfile(config{4}.statsavedir,'Waves_detections');

if ~isfolder(detect_stats_dir)
    mkdir(detect_stats_dir);
end

%store p-values in a single structure
Stats.pval_waves=p_val_waves;
Stats.adj_pval_waves=adj_pval_waves;

Stats.pval_trials=p_val_trials;
Stats.adj_pval_trials=adj_pval_trials;


save(fullfile(detect_stats_dir,'tests_mean_calculated_data.mat'),'Stats');


