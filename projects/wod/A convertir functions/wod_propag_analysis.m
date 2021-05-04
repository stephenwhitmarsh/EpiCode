function extracted_data= wod_propag_analysis(cfg,force)


detectiondatapath= fullfile(cfg{4}.datasavedir,'Detection');


fname_out=fullfile(detectiondatapath,'calculated_data.mat'); 
if exist(fname_out, 'file') && force == false
    load(fname_out, 'extracted_data');
    return
end


%% loading data
temp= load(fullfile(detectiondatapath,'stats_all.mat'));
stats_all=temp.stats_all;
clear temp

for irat= 1:size(cfg,2)
    
    if isempty(cfg{irat})
        continue
    end
    
    
    irat_name=sprintf('Rat_%i',irat);
    for iwave=["WoD" "WoR"]
    for itime= ["peak_time" "min_slope_time" "start_time"]
        for itrial= 1:size(stats_all.(iwave).peak_time,2)
            %% WOD Find origin and store origin depth and timings
            A=stats_all{irat}.(iwave).(itime)(:,itrial);
            B= stats_all{irat}.Depth(:,itrial);
            %find minimum timing
            origin_time=min(A);
            %store minimum timing in structure
            Data_waves.(iwave).origin_time.(itime)(irat,itrial)=origin_time;
            
            %find index of origin and get origin depth
            idx_origin=find(A==origin_time);
            origin_depth=stats_all{irat}.Depth(idx_origin,itrial);
            %store origin depth
            Data_waves.(iwave).origin_depth.(itime)(irat,itrial)=origin_depth;
            
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
                Data_waves.(iwave).speed.(itime).up(irat,itrial)=nanmean(Speed_up);
            end
            if bool_down==1
                 Data_waves(iwave).speed.(itime).down(irat,itrial)=nanmean(Speed_down);
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
            
            Speed_instant{irat}.(iwave).(itime)(:,itrial)=Speed_instan;
            
            
            
            clear A origin_time origin_depth bool_up bool_down Speed_up Speed_down Speed_instan
            
        end %itrial
    end %itime
    end %iwave
end %irat

%% Replace all 0 values by NaNs
for iwave=["WoD" "WoR"]
    for iana= ["origin_time" "origin_depth"]
        for itime= ["peak_time" "min_slope" "start_time"]
            idx=find(Data_waves{irat}.(iwave).(iana).(itime)(:,:)==0);
            Data_waves{irat}.(iwave).(iana).(itime)(idx)=NaN;
        end %itime
    end %iana
end %iwave

save(fullfile(detectiondatapath,'calculated_data.mat'),'Data_waves');
save(fullfile(detectiondatapath,'instant_speed.mat'),'Speed_instant')

%% Make statiscal analysis between trials

% origin time and depth
for iwave=["WOD" "WOR"]
    for iana= ["origin_time" "origin_depth"]
        count=0;
        for itime= ["peak_time" "min_slope" "start_time"]
            count= count+1;
            %separate trials as different arrays
            first_trial=Data_waves{irat}.(iwave).(iana).(itime)(:,1);
            second_trial=Data_waves{irat}.(iwave).(iana).(itime)(:,2);
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


