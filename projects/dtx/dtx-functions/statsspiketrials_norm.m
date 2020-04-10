function statsspiketrials_norm(cfg,SpikeTrials,ipart,ilabel, input_temp)

cfgtemp             = [];
cfgtemp.bins        = [0:0.005:150];
cfg.latency         = 'max';
cfg.keeptrials      = 'yes';
cfg.param           = 'coeffvar';
isidata = ft_spike_isi(cfgtemp,SpikeTrials{ipart}{ilabel});


%Ajouter CV2, CV, burst index, doublets de PA, baseline

% trialtime : debut, fin
% trial : index pour chaque spike du trial auquel il appartient
% time : time de chaque spike
% isidata.isi{itemplate} : isi entre chaque chaque spike. 1er spike du trial :ISI = NaN
smooth_cste = 0.1;

%TEMP, TO CHECK WHAT FREQLIM IS THE BEST
freqlim_list = {[0,10],[0,20],[0,30],[0,40],[0,50],[0,100],[0,200]}; %FIXME, choose the best
freqlim = freqlim_list(input_temp);
cfg.imagesavedir = fullfile(cfg.imagesavedir, ['freqlim_[',strrep(num2str(freqlim{1}), ' ', '_'),']']);
% check if images directory exists, if not create
if ~isfolder(cfg.imagesavedir)
    ft_notice('creating directory %s', cfg.imagesavedir);
    mkdir(cfg.imagesavedir);
end


for itemplate = 1:size(SpikeTrials{ipart}{ilabel}.label,2)
    
    %% One figure per template
    trialkept_index = [];
%     fig1 = figure(1);
%     hold;
%     fig2 = figure(2);
%     hold;
%     fig3 = figure(3);
%     hold;
%     fig4 = figure(4);
%     hold;
%     fig5 = figure(5);
%     hold;
    for itrial = 1:size(SpikeTrials{ipart}{ilabel}.trialinfo,1)
        spike_index = [];
        spike_index = (SpikeTrials{ipart}{ilabel}.trial{itemplate} == itrial); %first is NaN
        if sum(spike_index)>1 %if there is spikes in the trial
            isi_temp                                    = isidata.isi{itemplate}(spike_index);
            time_temp                                   = SpikeTrials{ipart}{ilabel}.time{itemplate}(spike_index); %time of the ISI in the trial
            isiTrials{itemplate}.trial{itrial}.isi            = isi_temp(2:end); %remove first NaN
            isiTrials{itemplate}.trial{itrial}.time           = time_temp(2:end); %remove first NaN
            %Normalisation : remove offset and divide by length
            isiTrials{itemplate}.trial{itrial}.time_norm      = (isiTrials{itemplate}.trial{itrial}.time - SpikeTrials{ipart}{ilabel}.trialtime(itrial,1)) ./ (SpikeTrials{ipart}{ilabel}.trialtime(itrial,2) - (SpikeTrials{ipart}{ilabel}.trialtime(itrial,1)));
            
            %smooth values
            for itime=1:1/smooth_cste %smooth_cste = time between 0 and 1 of the windows for smoothing ISI
                isi_smooth_index                            = [];
                isi_smooth_index                            = (isiTrials{itemplate}.trial{itrial}.time_norm > smooth_cste*(itime-1) & isiTrials{itemplate}.trial{itrial}.time_norm < smooth_cste*(itime));
                if sum(isi_smooth_index) >0
                    %Get ISI
                    isi_temp = [];
                    isi_temp = isiTrials{itemplate}.trial{itrial}.isi(isi_smooth_index);
                    isiTrials{itemplate}.trial{itrial}.timesmooth(itime)   =  smooth_cste*itime;
                    %Freq
                    isiTrials{itemplate}.trial{itrial}.isismooth(itime)    =  nanmean(isi_temp);
                    isiTrials{itemplate}.trial{itrial}.freqsmooth(itime)   =  1 / isiTrials{itemplate}.trial{itrial}.isismooth(itime);
                    %CV2, CV, FF
                    if sum(isi_smooth_index)>2
                        cv2_temp = [];
                        for i = 1:length(isi_temp)-1
                            cv2_temp(i) = 2*abs(isi_temp(i)-isi_temp(i+1))/(isi_temp(i)+isi_temp(i+1));
                        end
                        isiTrials{itemplate}.trial{itrial}.cv2(itime)      = nanmean(cv2_temp);
                        isiTrials{itemplate}.trial{itrial}.cv(itime)      = nanstd(isi_temp) / nanmean(isi_temp);
                        isiTrials{itemplate}.trial{itrial}.ff(itime)      = nanstd(isi_temp)^2 / nanmean(isi_temp);
                        isiTrials{itemplate}.trial{itrial}.burstindex(itime)       = sum(isi_temp<0.010)/sum(isi_temp<0.100);
                    else
                        isiTrials{itemplate}.trial{itrial}.cv2(itime)      = NaN;
                        isiTrials{itemplate}.trial{itrial}.cv(itime)       = NaN;
                        isiTrials{itemplate}.trial{itrial}.ff(itime)       = NaN;
                        isiTrials{itemplate}.trial{itrial}.burstindex(itime)       = NaN;
                    end
                    
                    
                else
                    isiTrials{itemplate}.trial{itrial}.isismooth(itime)    =  0;
                    isiTrials{itemplate}.trial{itrial}.freqsmooth(itime)   =  0;
                    isiTrials{itemplate}.trial{itrial}.timesmooth(itime)   =  smooth_cste*itime;
                    isiTrials{itemplate}.trial{itrial}.cv2(itime)      = NaN;
                    isiTrials{itemplate}.trial{itrial}.cv(itime)       = NaN;
                    isiTrials{itemplate}.trial{itrial}.ff(itime)       = NaN;
                    isiTrials{itemplate}.trial{itrial}.burstindex(itime)       = NaN;
                    %keep data for avg
                end
            end
            if max(isiTrials{itemplate}.trial{itrial}.freqsmooth) > freqlim{1}(1) &&  max(isiTrials{itemplate}.trial{itrial}.freqsmooth) < freqlim{1}(2)%remove bad trials
                trialkept_index = [trialkept_index, itrial];
%                 figure(1)
%                 plot(isiTrials{itemplate}.trial{itrial}.timesmooth, isiTrials{itemplate}.trial{itrial}.freqsmooth);
%                 figure(2)
%                 plot(isiTrials{itemplate}.trial{itrial}.timesmooth, isiTrials{itemplate}.trial{itrial}.cv2);
%                 figure(3)
%                 plot(isiTrials{itemplate}.trial{itrial}.timesmooth, isiTrials{itemplate}.trial{itrial}.cv);
%                 figure(4)
%                 plot(isiTrials{itemplate}.trial{itrial}.timesmooth, isiTrials{itemplate}.trial{itrial}.ff);
%                 figure(5)
%                 plot(isiTrials{itemplate}.trial{itrial}.timesmooth, isiTrials{itemplate}.trial{itrial}.burstindex);
            end
        else
            %Plot zero at each point if there are no AP in the trial
            %isiTrials{itemplate}.trial{itrial} = [];
            % plot(smooth_cste*(1:1/smooth_cste), zeros(1,fix(1/smooth_cste)));
        end
    end %itrial
    
    %get data for averaging
    for itime=1:1/smooth_cste
        for itrial = 1:size(SpikeTrials{ipart}{ilabel}.trialinfo,1)
            if ismember(itrial, trialkept_index)
                freqsmooth_temp{itime}{itrial} = isiTrials{itemplate}.trial{itrial}.freqsmooth(itime);
                cv2_template{itime}{itrial} = isiTrials{itemplate}.trial{itrial}.cv2(itime);
                cv_template{itime}{itrial} = isiTrials{itemplate}.trial{itrial}.cv(itime);
                ff_template{itime}{itrial} = isiTrials{itemplate}.trial{itrial}.ff(itime);
                burstindex_template{itime}{itrial} = isiTrials{itemplate}.trial{itrial}.burstindex(itime);
            end
        end
    end
    
    %plot avg
    for itime=1:1/smooth_cste
        isiTrials{itemplate}.avg.meanfreqsmooth(itime)      = nanmean([freqsmooth_temp{itime}{:}]);
        isiTrials{itemplate}.avg.medfreqsmooth(itime)       = nanmedian([freqsmooth_temp{itime}{:}]);
        isiTrials{itemplate}.avg.meancv2(itime)             = nanmean([cv2_template{itime}{:}]);
        isiTrials{itemplate}.avg.medcv2(itime)              = nanmedian([cv2_template{itime}{:}]);
        isiTrials{itemplate}.avg.meancv(itime)              = nanmean([cv_template{itime}{:}]);
        isiTrials{itemplate}.avg.medcv(itime)               = nanmedian([cv_template{itime}{:}]);
        isiTrials{itemplate}.avg.meanff(itime)              = nanmean([ff_template{itime}{:}]);
        isiTrials{itemplate}.avg.medff(itime)               = nanmedian([ff_template{itime}{:}]);
        isiTrials{itemplate}.avg.meanburstindex(itime)      = nanmean([burstindex_template{itime}{:}]);
        isiTrials{itemplate}.avg.medburstindex(itime)       = nanmedian([burstindex_template{itime}{:}]);
        %             isiTrials{itemplate}.avg.freqsmooth(itime)     = 1 / isiTrials{itemplate}.avg.isismooth(itime);
        isiTrials{itemplate}.avg.timesmooth(itime)          = smooth_cste*itime;
    end
    
%     figure(1)
%     plot(isiTrials{itemplate}.avg.timesmooth, isiTrials{itemplate}.avg.meanfreqsmooth, 'r', 'LineWidth', 2);
%     plot(isiTrials{itemplate}.avg.timesmooth, isiTrials{itemplate}.avg.medfreqsmooth, 'b', 'LineWidth', 2);
%     figure(2)
%     plot(isiTrials{itemplate}.avg.timesmooth, isiTrials{itemplate}.avg.meancv2, 'r', 'LineWidth', 2);
%     plot(isiTrials{itemplate}.avg.timesmooth, isiTrials{itemplate}.avg.medcv2, 'b', 'LineWidth', 2);
%     figure(3)
%     plot(isiTrials{itemplate}.avg.timesmooth, isiTrials{itemplate}.avg.meancv, 'r', 'LineWidth', 2);
%     plot(isiTrials{itemplate}.avg.timesmooth, isiTrials{itemplate}.avg.medcv, 'b', 'LineWidth', 2);
%     figure(4)
%     plot(isiTrials{itemplate}.avg.timesmooth, isiTrials{itemplate}.avg.meanff, 'r', 'LineWidth', 2);
%     plot(isiTrials{itemplate}.avg.timesmooth, isiTrials{itemplate}.avg.medff, 'b', 'LineWidth', 2);
%     figure(5)
%     plot(isiTrials{itemplate}.avg.timesmooth, isiTrials{itemplate}.avg.meanburstindex, 'r', 'LineWidth', 2);
%     plot(isiTrials{itemplate}.avg.timesmooth, isiTrials{itemplate}.avg.medburstindex, 'b', 'LineWidth', 2);
%     
%     
%     set(fig1,'PaperOrientation','landscape');
%     set(fig1,'PaperUnits','normalized');
%     set(fig1,'PaperPosition', [0 0 1 1]);
%     print(fig1, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.name{ilabel},'_',SpikeTrials{ipart}{ilabel}.label{itemplate},'_FreqOverTrials.pdf']));
%     print(fig1, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.name{ilabel},'_',SpikeTrials{ipart}{ilabel}.label{itemplate},'_FreqOverTrials.png']),'-r600');
%     
%     set(fig2,'PaperOrientation','landscape');
%     set(fig2,'PaperUnits','normalized');
%     set(fig2,'PaperPosition', [0 0 1 1]);
%     print(fig2, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.name{ilabel},'_',SpikeTrials{ipart}{ilabel}.label{itemplate},'_CV2OverTrials.pdf']));
%     print(fig2, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.name{ilabel},'_',SpikeTrials{ipart}{ilabel}.label{itemplate},'_CV2OverTrials.png']),'-r600');
% 
%     set(fig3,'PaperOrientation','landscape');
%     set(fig3,'PaperUnits','normalized');
%     set(fig3,'PaperPosition', [0 0 1 1]);
%     print(fig3, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.name{ilabel},'_',SpikeTrials{ipart}{ilabel}.label{itemplate},'_CVOverTrials.pdf']));
%     print(fig3, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.name{ilabel},'_',SpikeTrials{ipart}{ilabel}.label{itemplate},'_CVOverTrials.png']),'-r600');
% 
%     set(fig4,'PaperOrientation','landscape');
%     set(fig4,'PaperUnits','normalized');
%     set(fig4,'PaperPosition', [0 0 1 1]);
%     print(fig4, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.name{ilabel},'_',SpikeTrials{ipart}{ilabel}.label{itemplate},'_FanoFactorOverTrials.pdf']));
%     print(fig4, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.name{ilabel},'_',SpikeTrials{ipart}{ilabel}.label{itemplate},'_FanoFactorOverTrials.png']),'-r600');
% 
%     set(fig5,'PaperOrientation','landscape');
%     set(fig5,'PaperUnits','normalized');
%     set(fig5,'PaperPosition', [0 0 1 1]);
%     print(fig5, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.name{ilabel},'_',SpikeTrials{ipart}{ilabel}.label{itemplate},'_BurstIndexOverTrials.pdf']));
%     print(fig5, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.name{ilabel},'_',SpikeTrials{ipart}{ilabel}.label{itemplate},'_BurstIndexOverTrials.png']),'-r600');

    close all
end %itemplate
%%
% get data for pooling over templates
for itime=1:1/smooth_cste
    for itemplate = 1:1:size(SpikeTrials{ipart}{ilabel}.label,2)
        freqovertrials_mean{itemplate}(itime)          = isiTrials{itemplate}.avg.meanfreqsmooth(itime);
        freqovertrials_med{itemplate}(itime)           = isiTrials{itemplate}.avg.medfreqsmooth(itime);        
        cv2overtrials_mean{itemplate}(itime)           = isiTrials{itemplate}.avg.meancv2(itime);
        cv2overtrials_med{itemplate}(itime)            = isiTrials{itemplate}.avg.medcv2(itime);        
        cvovertrials_mean{itemplate}(itime)            = isiTrials{itemplate}.avg.meancv(itime);
        cvovertrials_med{itemplate}(itime)             = isiTrials{itemplate}.avg.medcv(itime);
        ffovertrials_mean{itemplate}(itime)            = isiTrials{itemplate}.avg.meanff(itime);
        ffovertrials_med{itemplate}(itime)             = isiTrials{itemplate}.avg.medff(itime);
        burstindexovertrials_mean{itemplate}(itime)    = isiTrials{itemplate}.avg.meanburstindex(itime);
        burstindexovertrials_med{itemplate}(itime)     = isiTrials{itemplate}.avg.medburstindex(itime);
    end
end

%Freq
fig = figure;
hold
for  itemplate = 1:size(SpikeTrials{ipart}{ilabel}.label,2)
    plot(smooth_cste*(1:1/smooth_cste),freqovertrials_mean{itemplate});
end

set(fig,'PaperOrientation','landscape');
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition', [0 0 1 1]);
print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.name{ilabel},'_','Alllates','_FreqOverTrials_Mean.pdf']));
print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.name{ilabel},'_','All_Templates','_FreqOverTrials_Mean.png']),'-r600');


fig = figure;
hold
for  itemplate = 1:size(SpikeTrials{ipart}{ilabel}.label,2)
    plot(smooth_cste*(1:1/smooth_cste),freqovertrials_med{itemplate});
end

set(fig,'PaperOrientation','landscape');
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition', [0 0 1 1]);
print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.name{ilabel},'_','All_Templates','_FreqOverTrials_Med.pdf']));
print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.name{ilabel},'_','All_Templates','_FreqOverTrials_Med.png']),'-r600');

fig = figure;
hold
for  itemplate = 1:size(SpikeTrials{ipart}{ilabel}.label,2)
    plot(smooth_cste*(1:1/smooth_cste),freqovertrials_med{itemplate});
end

ylim([0 5]);
set(fig,'PaperOrientation','landscape');
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition', [0 0 1 1]);
print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.name{ilabel},'_','All_Templates','_FreqOverTrials_Med_littlescale.pdf']));
print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.name{ilabel},'_','All_Templates','_FreqOverTrials_Med_littlescale.png']),'-r600');

%CV2
fig = figure;
hold
for  itemplate = 1:size(SpikeTrials{ipart}{ilabel}.label,2)
    plot(smooth_cste*(1:1/smooth_cste),cv2overtrials_mean{itemplate});
end

set(fig,'PaperOrientation','landscape');
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition', [0 0 1 1]);
print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.name{ilabel},'_','All_Templates','_CV2OverTrials_Mean.pdf']));
print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.name{ilabel},'_','All_Templates','_CV2OverTrials_Mean.png']),'-r600');


fig = figure;
hold
for  itemplate = 1:size(SpikeTrials{ipart}{ilabel}.label,2)
    plot(smooth_cste*(1:1/smooth_cste),cv2overtrials_med{itemplate});
end

set(fig,'PaperOrientation','landscape');
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition', [0 0 1 1]);
print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.name{ilabel},'_','All_Templates','_CV2OverTrials_Med.pdf']));
print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.name{ilabel},'_','All_Templates','_CV2OverTrials_Med.png']),'-r600');

%CV
fig = figure;
hold
for  itemplate = 1:size(SpikeTrials{ipart}{ilabel}.label,2)
    plot(smooth_cste*(1:1/smooth_cste),cvovertrials_mean{itemplate});
end

set(fig,'PaperOrientation','landscape');
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition', [0 0 1 1]);
print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.name{ilabel},'_','All_Templates','_CVOverTrials_Mean.pdf']));
print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.name{ilabel},'_','All_Templates','_CVOverTrials_Mean.png']),'-r600');


fig = figure;
hold
for  itemplate = 1:size(SpikeTrials{ipart}{ilabel}.label,2)
    plot(smooth_cste*(1:1/smooth_cste),cvovertrials_med{itemplate});
end

set(fig,'PaperOrientation','landscape');
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition', [0 0 1 1]);
print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.name{ilabel},'_','All_Templates','_CVOverTrials_Med.pdf']));
print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.name{ilabel},'_','All_Templates','_CVOverTrials_Med.png']),'-r600');

% FF
fig = figure;
hold
for  itemplate = 1:size(SpikeTrials{ipart}{ilabel}.label,2)
    plot(smooth_cste*(1:1/smooth_cste),ffovertrials_mean{itemplate});
end

set(fig,'PaperOrientation','landscape');
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition', [0 0 1 1]);
print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.name{ilabel},'_','All_Templates','_FanoFactorOverTrials_Mean.pdf']));
print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.name{ilabel},'_','All_Templates','_FanoFactorOverTrials_Mean.png']),'-r600');


fig = figure;
hold
for  itemplate = 1:size(SpikeTrials{ipart}{ilabel}.label,2)
    plot(smooth_cste*(1:1/smooth_cste),ffovertrials_med{itemplate});
end

set(fig,'PaperOrientation','landscape');
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition', [0 0 1 1]);
print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.name{ilabel},'_','All_Templates','_FanoFactorOverTrials_Med.pdf']));
print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.name{ilabel},'_','All_Templates','_FanoFactorOverTrials_Med.png']),'-r600');

%Burst Index
fig = figure;
hold
for  itemplate = 1:size(SpikeTrials{ipart}{ilabel}.label,2)
    plot(smooth_cste*(1:1/smooth_cste),burstindexovertrials_mean{itemplate});
end

set(fig,'PaperOrientation','landscape');
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition', [0 0 1 1]);
print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.name{ilabel},'_','All_Templates','_burstindexOverTrials_Mean.pdf']));
print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.name{ilabel},'_','All_Templates','_burstindexOverTrials_Mean.png']),'-r600');


fig = figure;
hold
for  itemplate = 1:size(SpikeTrials{ipart}{ilabel}.label,2)
    plot(smooth_cste*(1:1/smooth_cste),burstindexovertrials_med{itemplate});
end

set(fig,'PaperOrientation','landscape');
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition', [0 0 1 1]);
print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.name{ilabel},'_','All_Templates','_burstindexOverTrials_Med.pdf']));
print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.name{ilabel},'_','All_Templates','_burstindexOverTrials_Med.png']),'-r600');

close all



end

