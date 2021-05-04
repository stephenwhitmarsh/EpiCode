function freq_data= wod_tfr_extractdata(cfg,force)

fname_out = fullfile(cfg{4}.datasavedir,'freq_data', sprintf('freq_data_allrat.mat'));
if exist(fname_out, 'file') && force == false
    load(fname_out, 'stats_all');
    return
end




freqstatpath= fullfile(cfg{4}.statsavedir,'freq_data');
analysis_names={'timefreq_wod','timefreq_wod_timenorm','timefreq_wod_blcorrected','timefreq_wod_timenorm_blcorrected'};

%% Prepare structure to analyze

for idata= 1: size(analysis_names,2)
    for irat= 4:size(cfg,2)
        cfgcommon = cfg{4}; %same common config for all rats
        
        fprintf('Load timefreq data for rat %d/%d\n', irat,size(cfg,2));
        data_temp = load(fullfile(cfg{irat}.datasavedir,[cfg{irat}.prefix,analysis_names{idata},'.mat']));
        
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
                    cfgtemp.frequency= cfg{irat}.timefreq.(iband);
                    %cfgtemp.avgoverfreq = 'yes';
                    pow.(iband).(chan_name)= ft_selectdata(cfgtemp, data_rat.(chan_name));
                    %average over frequency
                    pow.(iband).(chan_name).powspctrm = permute(pow.(iband).(chan_name).powspctrm, [2 3 1]); %put first dimension as last so it is deleted
                    pow.(iband).(chan_name).powspctrm = nanmean(pow.(iband).(chan_name).powspctrm,1);
                    pow.(iband).(chan_name).time= pow.(iband).(chan_name).time;
                    
                    if iband == "LF"
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
                    
                    freq_data{irat}.(analysis_names{idata}).peak_time.(iband)(ichan,itrial)=  t_peak.(iband).(chan_name);
                    freq_data{irat}.(analysis_names{idata}).peak_value.(iband)(ichan,itrial)=  v_peak.(iband).(chan_name);

                    
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
                fname= fullfile(cfg{4}.imagesavedir_data{2},'detection',sprintf('%sWOD%i_%s_%s',cfg{irat}.prefix,itrial,iband,analysis_names{idata}));
                dtx_savefigure(fig,fname,'png','pdf','close');
            end %iband
        end %itrial
    end %irat
end %idata

%% Save data

if ~isfolder(fullfile(cfg{4}.datasavedir,'freq_data'))
    mkdir(fullfile(cfg{4}.datasavedir,'freq_data'));
end

save(fname_out,'freq_data');

end %function