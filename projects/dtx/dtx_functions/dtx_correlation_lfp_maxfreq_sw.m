%stats = obteus avec dtx_project_spike, modifiées avec doc excl

%plot freq en fonction de max LFP
%ou plot temps max freq en fonction de temps max LFP

ipart = 1;
ievent = 1;
rat_list = 1:5;

% Load LFP for each rat
for irat = rat_list
    [MuseStruct]                    = readMuseMarkers(config{irat}, false);
    [MuseStruct]                    = alignMuseMarkers(config{irat},MuseStruct, false);
    LFP{irat}                       = readLFP(config{irat},MuseStruct,false);
    LFP{irat}                       = removetrials_MuseMarkers(config{irat}, LFP{irat}, MuseStruct);
end

%Compute amplitude LFP for each channel, for each trial
%select only the peaks > 3*SD of BL
for irat = rat_list
    %lp filter data 
    cfgtemp = [];
    cfgtemp.lpfilter = 'yes';
    cfgtemp.lpfreq = 10;
    data = ft_preprocessing(cfgtemp,LFP{irat}{ipart}{ievent});
    for ichan = 1:size(LFP{irat}{ipart}{ievent}.label,1)
        ft_progress('init', 'text',     sprintf('Rat %d/%d, Chan %d/%d:', irat,5,ichan, size(data.label,1)));
        for itrial = 1:size(LFP{irat}{ipart}{ievent}.trial,2)
            
            ft_progress(itrial/size(data.trial,2), 'computing amplitude for trial %d from %d', itrial, size(data.trial,2));
            amplitude_chanlabel{irat}{ichan} = data.label{ichan};
            
            %compute bl
            bl_idx = data.time{itrial} >= -2 & data.time{itrial} <= -1;
            bl     = mean(data.trial{itrial}(ichan,bl_idx));   
            std_bl  = std(data.trial{itrial}(ichan,bl_idx));   
            
            %compute peak
            ac_idx = data.time{itrial} >= -0.5 & data.time{itrial} <= 1;
            [peak, loc] = findpeaks(data.trial{itrial}(ichan,ac_idx),data.time{itrial}(ac_idx),'NPeaks',1,'SortStr','descend','WidthReference','Halfheight'); %Npeaks : max nr of peaks/ SortStr : peak sorting : descend = from largest to smallest

            if isempty(peak)
                amplitude{irat}{ichan}(itrial)   = NaN;
                amplitude_loc{irat}{ichan}(itrial) = NaN;
                continue
            end
            if peak < 3*std_bl
                amplitude{irat}{ichan}(itrial)   = NaN;
                amplitude_loc{irat}{ichan}(itrial) = NaN;
                continue
            end
            
            amplitude{irat}{ichan}(itrial)   = peak-bl;
            amplitude_loc{irat}{ichan}(itrial) = loc;
        end 
        ft_progress('close');
    end
end

% get max sdf value for each trial, for each neuron
for irat = rat_list
    sdf = stats{irat}{ipart}.SlowWave.sdfdata;
    %             figure;hold;
    for i_unit = 1:size(stats{irat}{ipart}.SlowWave.label,2)
        spike_maxchan{irat}(i_unit) = stats{irat}{ipart}.Interictal.maxchan{i_unit};
        for itrial = 1:size(sdf.trial,2)
            time_sel = sdf.time{itrial} >= -0.5 & sdf.time{itrial}<= 1;
            [peak loc] = findpeaks(sdf.trial{itrial}(i_unit,time_sel),sdf.time{itrial}(time_sel),'NPeaks',1,'SortStr','descend'); %Npeaks : max nr of peaks/ SortStr : peak sorting : descend = from largest to smallest
            if ~isempty(peak)
                maxfreq{irat}{i_unit}(itrial) = peak;
                maxfreq_loc{irat}{i_unit}(itrial) = loc;
            else
                maxfreq{irat}{i_unit}(itrial) = NaN;
                maxfreq_loc{irat}{i_unit}(itrial) = NaN;
            end
            %                 plot(sdf.trial{itrial}(i_unit,time_sel), 'k');
            %                 scatter(loc(i_unit), peak(i_unit), 'xr');
        end
    end
end

%get lfp amplitude of the chan for each neuron
% and sua.mua, pn/in
for irat = rat_list
    for i_unit = 1:size(stats{irat}{ipart}.SlowWave.label,2)
        if contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}, 'noise')
            alldata.swamplitude{irat}{i_unit}     = NaN;
            alldata.group{irat}{i_unit}           = 'noise';
            continue
        end
        alldata.label{irat}{i_unit}           = stats{irat}{ipart}.(config{irat}.spike.baselinename).label{i_unit};
        alldata.group{irat}{i_unit}           = stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit};
        alldata.celltype{irat}{i_unit}        = stats{irat}{ipart}.(config{irat}.spike.baselinename).celltype{i_unit};
        alldata.code_slowwave{irat}{i_unit}   = stats{irat}{ipart}.SlowWave.code_slowwave{i_unit};
        alldata.maxchan{irat}{i_unit}         = sprintf('E%.2dLFP',stats{irat}{ipart}.(config{irat}.spike.baselinename).maxchan{i_unit});
        chansel = find(strcmp(amplitude_chanlabel{irat}, alldata.maxchan{irat}{i_unit}));
        m1g_sel = find(strcmp(amplitude_chanlabel{irat}, 'ECoGM1G'));
        if irat == 2, m1g_sel = find(strcmp(amplitude_chanlabel{irat}, 'ECoGS1')); end %erreur pendant l'acquisition
        alldata.swamplitude{irat}{i_unit}     = amplitude{irat}{chansel};
        alldata.maxfreq{irat}{i_unit}         = maxfreq{irat}{i_unit};        
        alldata.swamplitude_loc{irat}{i_unit} = amplitude_loc{irat}{chansel};
        alldata.maxfreq_loc{irat}{i_unit}     = maxfreq_loc{irat}{i_unit};
        alldata.m1g_peak_loc{irat}{i_unit}    = amplitude_loc{irat}{m1g_sel};
    end
end
        
%% plot maxfreq en fonction de amplitude
% x = lfp amplitude, y = maxfreq
% one color per unit
% plot avec only sua, et only mua, et les deux
% plot avec only pn, et only in
% plot en fonction du code slowwave aussi

% compute linear fit, rho and pval

figure;hold;setfig();
count_pos = 0;
count_neg = 0;
C = linspecer(50,'qualitative');
i=0;
for irat = rat_list
    for i_unit = 1:size(alldata.label{irat},2)
        if contains(alldata.group{irat}{i_unit}, 'noise')
            continue
        end
        if contains(alldata.group{irat}{i_unit}, 'mua')
            dofill = false;
%             continue
        end
        if ~contains(alldata.group{irat}{i_unit}, 'mua') %sua
            dofill = true;
%             continue
        end
        if contains(alldata.celltype{irat}{i_unit}, 'pn')
            plottype = '^k';
%             continue
        end
        if contains(alldata.celltype{irat}{i_unit}, 'in')
            plottype = 'ok';
%             continue
        end
        if isempty(alldata.maxfreq{irat}{i_unit}) || isempty(alldata.swamplitude{irat}{i_unit})
            continue
        end
%         i=i+1;
%         color=C(i,:);
        x = alldata.swamplitude{irat}{i_unit};
        y = alldata.maxfreq{irat}{i_unit};
%         if dofill
%             scatter(x, y,plottype, 'filled','MarkerEdgeColor',color,'MarkerFaceColor',color);
%         else
%             scatter(x, y,plottype,'MarkerEdgeColor',color);
%         end
%         
        sel = ~isnan(y);
        x=x(sel);
        y=y(sel);
        
        coeffs_fit = polyfit(x,y,1);
        xFit{irat}{i_unit} = linspace(min(x), max(x), 500);
        yFit{irat}{i_unit} = polyval(coeffs_fit, xFit{irat}{i_unit});
        yFit_norm{irat}{i_unit} = normalize(yFit{irat}{i_unit});
%         plot(xFit{irat}{i_unit}, yFit{irat}{i_unit},'b','LineWidth',2);
        
        [rho{irat}{i_unit},pval{irat}{i_unit}] = corr(x',y','Type','Spearman');
        
        if pval{irat}{i_unit} < 0.001
            count_pos = count_pos+1;
            color = 'b';
        elseif pval{irat}{i_unit} < 0.05
            count_pos = count_pos+1;
            color = 'r';
        else
            count_neg = count_neg+1;
            color = [0.6 0.6 0.6];
        end
        plot(xFit{irat}{i_unit}, yFit{irat}{i_unit},'Color',color,'LineWidth',2);
    end
end
% xlabel('SW amplitude (µV)');
% ylabel('Max freq during SW (µV)');

ylim([0 150]);
setfig();

%plot rho
figure;hold;setfig();
for irat = rat_list
    for i_unit = 1:size(alldata.label{irat},2)
        if contains(alldata.group{irat}{i_unit}, 'noise')
            continue
        end
        if contains(alldata.group{irat}{i_unit}, 'mua')
            dofill = false;
%             continue
        end
        if ~contains(alldata.group{irat}{i_unit}, 'mua') %sua
            dofill = true;
%             continue
        end
        if contains(alldata.celltype{irat}{i_unit}, 'pn')
            plottype = '^k';
            %             continue
        end
        if contains(alldata.celltype{irat}{i_unit}, 'in')
            plottype = 'ok';
            %             continue
        end
%         y = rho{irat}{i_unit};
        y = pval{irat}{i_unit};
        if dofill
            scatter(rand, y,plottype, 'filled');
        else
            scatter(rand, y,plottype);
        end
    end
end
xlim([-1 2]);
set(gca,'Yscale', 'linear');
plot([-1 2], [0.05 0.05], '--r');
%plot pval


% clear xFit yFit rho pval

%% plot délais maxfreq en fonction de délai LFP, ou EEG M1G
figure;hold;
% C = linspecer(50,'qualitative');
i=0;
for irat = rat_list
    for i_unit = 1:size(alldata.label{irat},2)
        if contains(alldata.group{irat}{i_unit}, 'noise')
            continue
        end
        if contains(alldata.group{irat}{i_unit}, 'mua')
            dofill = false;
%             continue
        end
        if ~contains(alldata.group{irat}{i_unit}, 'mua') %sua
            dofill = true;
%             continue
        end
        if contains(alldata.celltype{irat}{i_unit}, 'pn')
            plottype = '^k';
%             continue
        end
        if contains(alldata.celltype{irat}{i_unit}, 'in')
            plottype = 'ok';
%             continue
        end
        if isempty(alldata.maxfreq{irat}{i_unit}) || isempty(alldata.swamplitude{irat}{i_unit})
            continue
        end
        i=i+1;
        color=C(i,:);
%         x = rand(size(alldata.maxfreq_loc{irat}{i_unit}));
        x = rand;
%         y = alldata.maxfreq_loc{irat}{i_unit} - alldata.swamplitude_loc{irat}{i_unit};
        y = nanmean(alldata.maxfreq_loc{irat}{i_unit}) - nanmean(alldata.m1g_peak_loc{irat}{i_unit});
%         if dofill
%             scatter(x, y,plottype, 'filled');%,'MarkerEdgeColor',color,'MarkerFaceColor',color);
%         else
%             scatter(x, y,plottype);%,'MarkerEdgeColor',color);
%         end
        if dofill
            scatter(x, y,plottype, 'filled');%,'MarkerEdgeColor',color,'MarkerFaceColor',color);
        else
            scatter(x, y,plottype);%,'MarkerEdgeColor',color);
        end
    end
end
% xlabel('SW peak time (µV)');
ylabel('Distribution of maxfreq and sw delays');

