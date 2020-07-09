ipart = 1;
ievent = 1;

% Load LFP for each rat
for irat = 1:5
    [MuseStruct]                    = readMuseMarkers(config{irat}, false);
    [MuseStruct]                    = alignMuseMarkers(config{irat},MuseStruct, false);
    LFP{irat}                       = readLFP(config{irat},MuseStruct,false);
    LFP{irat}                       = removetrials_MuseMarkers(config{irat}, LFP{irat}, MuseStruct);
end

%Compute amplitude LFP for each channel, for each trial
for irat = 1:5
    %lp filter data 
    cfgtemp = [];
    cfgtemp.lpfilter = 'yes';
    cfgtemp.lpfreq = 10;
    data = ft_preprocessing(cfgtemp,LFP{irat}{ipart}{ievent});
    for ichan = 1:size(LFP{irat}{ipart}{ievent}.label,1)
        ft_progress('init', 'text',     sprintf('Rat %d/%d, Chan %d/%d:', irat,5,ichan, size(data.label,1)));
        for itrial = 1:size(LFP{irat}{ipart}{ievent}.trial,2)
            
            ft_progress(itrial/size(data.trial,2), 'computing amplitude for trial %d from %d', itrial, size(data.trial,2));
            
            %compute bl
            bl_idx = data.time{itrial} >= -2 & data.time{itrial} <= -1;
            bl     = mean(data.trial{itrial}(ichan,bl_idx));       
            
            %compute peak
            ac_idx = data.time{itrial} >= -0.5 & data.time{itrial} <= 1;
            peak = findpeaks(data.trial{itrial}(ichan,ac_idx),'NPeaks',1,'SortStr','descend','WidthReference','Halfheight'); %Npeaks : max nr of peaks/ SortStr : peak sorting : descend = from largest to smallest

            if ~isempty(peak)
                amplitude{irat}{ichan}(itrial)   = peak-bl;
            else
                amplitude{irat}{ichan}(itrial)   = NaN;
            end
            amplitude_chanlabel{irat}{ichan} = data.label{ichan};
        end 
        ft_progress('close');
    end
end

% get max sdf value for each trial, for each neuron
for irat = 1:5
    sdf = stats{irat}{ipart}.SlowWave.sdfdata;
    %             figure;hold;
    for i_unit = 1:size(stats{irat}{ipart}.SlowWave.label,2)
        spike_maxchan{irat}(i_unit) = stats{irat}{ipart}.Interictal.maxchan{i_unit};
        for itrial = 1:size(sdf.trial,2)
            time_sel = sdf.time{itrial} >= -0.5 & sdf.time{itrial}<= 1;
            peak = findpeaks(sdf.trial{itrial}(i_unit,time_sel),'NPeaks',1,'SortStr','descend'); %Npeaks : max nr of peaks/ SortStr : peak sorting : descend = from largest to smallest
            if ~isempty(peak)
                maxfreq{irat}{i_unit}(itrial) = peak;
            else
                maxfreq{irat}{i_unit}(itrial) = NaN;
            end
            %                 plot(sdf.trial{itrial}(i_unit,time_sel), 'k');
            %                 scatter(loc(i_unit), peak(i_unit), 'xr');
        end
    end
end

%get lfp amplitude of the chan for each neuron
% and sua.mua, pn/in
for irat = 1:5
    for i_unit = 1:size(stats{irat}{ipart}.SlowWave.label,2)
        if contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}, 'noise')
            continue
        end
        alldata.label{irat}{i_unit}           = stats{irat}{ipart}.(config{irat}.spike.baselinename).label{i_unit};
        alldata.group{irat}{i_unit}           = stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit};
        alldata.celltype{irat}{i_unit}        = stats{irat}{ipart}.(config{irat}.spike.baselinename).celltype{i_unit};
        alldata.maxchan{irat}{i_unit}         = sprintf('E%.2dLFP',stats{irat}{ipart}.(config{irat}.spike.baselinename).maxchan{i_unit});
        chansel = find(strcmp(amplitude_chanlabel{irat}, alldata.maxchan{irat}{i_unit}));
        alldata.swamplitude{irat}{i_unit}     = amplitude{irat}{chansel};
        alldata.maxfreq{irat}{i_unit}         = maxfreq{irat}{i_unit};
    end
end
        
%plot 
% x = lfp amplitude, y = maxfreq
% one color per unit
% plot avec only sua, et only mua, et les deux
% plot avec only pn, et only in
figure;hold;
for irat = 1:5
    for i_unit = 1:size(alldata.label{irat},2)
        if contains(alldata.group{i_unit}, 'noise')
            continue
        end
        if contains(alldata.group{i_unit}, 'mua')
            dofill = false;
%             continue
        end
        if ~contains(alldata.group{i_unit}, 'mua') %sua
            dofill = true;
%             continue
        end
        if contains(alldata.celltype{i_unit}, 'pn')
            plottype = '^k';
%             continue
        end
        if contains(alldata.celltype{i_unit}, 'in')
            plottype = 'ok';
%             continue
        end
        if dofill
            scatter(alldata.swamplitude{irat}{i_unit}, alldata.maxfreq{irat}{i_unit},plottype, 'filled');
        else
            scatter(alldata.swamplitude{irat}{i_unit}, alldata.maxfreq{irat}{i_unit},plottype);
        end
    end
end
