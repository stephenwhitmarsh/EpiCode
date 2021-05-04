function ordered_data=wod_fusion_data(data,cfg,force)

%load data
fname_out = fullfile(cfg{4}.datasavedir,'Detection', sprintf('wod_wavedetection_allprobes.mat'));
if exist(fname_out, 'file') && force == false
    load(fname_out, 'stats');
    return
end


%concat data
temp=load(fullfile(cfg{4}.datasavedir,'Detection', sprintf('wod_wavedetection_allrats.mat')));
stats=temp.stats;
clear temp


%% Arrange structure for 32 chan 

%FIXME: need WOD_data, WOR_data, Instan_Speed, Peak_time_freqband,
%peak_value_freqband


depth_start = 0;
depth_end = 2500; %µm
depth_step = 100;

depth_allrats = [];
wodtimings_all = [];
icol = 0;
for idepth = depth_start:depth_step:depth_end %step de 250 car les électrodes sont espacées de 250µm
    icol = icol+1;
    irow = 0;
    for ratname = string(fieldnames(Electrode_depth))'
        for itrial = 1:size(Electrode_depth.(ratname), 2)
            irow = irow+1;
            sel = abs(Electrode_depth.(ratname)(:, itrial) - idepth) < depth_step/2;
            if sum(sel) == 1
                depth_allrats(irow,icol) = Electrode_depth.(ratname)(sel, itrial);
                wodtimings_all(irow, icol) = WOD_data.timings.(ratname).peak(sel,itrial);
            elseif sum(sel) == 0
                depth_allrats(irow,icol) = nan;
                wodtimings_all(irow, icol) = nan;
            elseif sum(sel) > 1
                error('it should have only one electrode for one deepness');
            end
        end
    end
end

%interpolate missing data
for irat = 1:size(wodtimings_all, 1)
    startinterp = find(~isnan(wodtimings_all(irat,:)), 1, 'first');
    endinterp = find(~isnan(wodtimings_all(irat,:)), 1, 'last');
    wodtimings_all(irat,startinterp:endinterp) = fillmissing(wodtimings_all(irat,startinterp:endinterp), 'linear'); %ou pchip, ou spline
end
    
avg = mean(wodtimings_all, 1, 'omitnan');
std = std(wodtimings_all, 0, 1, 'omitnan');
figure; hold on;
% plot(avg,depth_start:depth_step:depth_end);
plot(depth_start:depth_step:depth_end, avg);
patch_std(depth_start:depth_step:depth_end, avg,std);
% patch_std(avg,depth_start:depth_step:depth_end, std);

figure; hold on
for irat = 1:size(wodtimings_all, 1)
    plot(wodtimings_all(irat, :) - nanmin(wodtimings_all(irat, :)), depth_start:depth_step:depth_end, 'k');
end

figure; hold on
for irat = 1:size(wodtimings_all, 1)
    toplot(irat , :) = wodtimings_all(irat, :) - nanmin(wodtimings_all(irat, :));
end


%% 

%load data for 16 and 32 chans

stats_concat.ratname = [];
depth_start = 0;
depth_end = 3200; %µm
depth_step = 100;

icol = 0;
for idepth = depth_start:depth_step:depth_end %step de 250 car les électrodes sont espacées de 250µm
    icol = icol+1;
    irow = 0;
    for irat = 1:size(stats_all, 2)
        for itrial = 1:size(stats_all{irat}.wod_peak_time, 2)
            irow = irow+1;
            sel = abs(stats_all{irat}.electrode_depth(:, itrial) - idepth) < depth_step/2;
            for ifield = string(fieldnames(stats_all{irat}))'
                if strcmp(ifield, 'ratname')
                    stats_concat.ratname = [stats_concat.ratname; string(stats_all{irat}.ratname)];
                    continue
                end
                if sum(sel) == 1
                    stats_concat.(ifield)(irow, icol) = stats_all{irat}.(ifield)(sel, itrial);
                    
                elseif sum(sel) == 0
                    stats_concat.(ifield)(irow, icol) = nan;
                elseif sum(sel) > 1
                    error('it should have only one electrode for one deepness');
                end
            end
        end
    end
end
%interpolate missing data
for irat = 1:size(wodtimings_all, 1)
    startinterp = find(~isnan(wodtimings_all(irat,:)), 1, 'first');
    endinterp = find(~isnan(wodtimings_all(irat,:)), 1, 'last');
    wodtimings_all(irat,startinterp:endinterp) = fillmissing(wodtimings_all(irat,startinterp:endinterp), 'linear'); %ou pchip, ou spline
end