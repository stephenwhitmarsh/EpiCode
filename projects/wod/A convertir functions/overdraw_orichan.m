if ispc
    addpath (genpath('\\lexport\iss01.charpier\analyses\wod\EpiCode\shared'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\wod\EpiCode\external'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\wod\EpiCode\templates'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\wod\EpiCode\projects\wod'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\wod\EpiCode\projects\dtx'))
    addpath \\lexport\iss01.charpier\analyses\wod\fieldtrip-20200607
    
elseif isunix
    addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/EpiCode/shared'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/EpiCode/external'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/EpiCode/templates'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/EpiCode/projects/wod'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/EpiCode/projects/dtx'))
    addpath /network/lustre/iss01/charpier/analyses/wod/fieldtrip-20200607
end

ft_defaults

config = wod_setparams;


%% Load data for Wod and WoR


% load each rat and plot 

for irat= 4: size(config,2)
    
    fprintf('Load LFP data for rat %d/%d\n', irat,size(config,2));
    temp=load(fullfile(config{irat}.datasavedir,[config{irat}.prefix,'LFP_WoD','.mat']));
    LFP = temp.LFP{1}.WoD;
    clear temp
    
    
    fprintf('Load Muse Markers for rat %d/%d\n', irat, size(config,2));
    temp= load(fullfile(config{irat}.datasavedir,[config{irat}.prefix,'MuseStruct','.mat']));
    MuseStruct =temp.MuseStruct{1};
    clear temp
    
%rename chans according to their real deepness.
        %the name is in cfg.LFP.channel, and it is renamed with the name at
        %the same index in cfg.LFP.rename
        %16 is surface, 1 is the deepest. 0 is the respi.
        n_chans = size(config{irat}.LFP.allchannel,2);
        for ichan = 1:n_chans
            if any(strcmp(config{irat}.LFP.channel,config{irat}.LFP.allchannel{ichan}))
                %search channel into config
                chan_idx = strcmp(config{irat}.LFP.channel,config{irat}.LFP.allchannel{ichan});
                new_name = config{irat}.LFP.rename{chan_idx};
                %search channel into LFP data to remane it
                chan_idx = strcmp(LFP.label, config{irat}.LFP.allchannel{ichan});
                LFP.label{chan_idx} = new_name;
            end
        end

%% Detect WoD for each trial

    %remove breathing and ekg channel
        cfgtemp         = [];
        cfgtemp.channel = {'all', '-E0', '-Respi', '-ECG','-Puff'};
        LFP             = ft_selectdata(cfgtemp, LFP);
        LFP_cleaned     = LFP; %save for later removing of artefacts
    
   % Separate LFP for trials
    for itrial= 1: size(LFP.trial,2)
         %select one trial
            cfgtemp         = [];
            cfgtemp.trials  = itrial;
            LFP_trial       = ft_selectdata(cfgtemp, LFP_cleaned);
            
            %filter lfp to better recognize WOD/WOR peak
            cfgtemp             = [];
            cfgtemp.lpfilter    = 'yes';
            cfgtemp.lpfilttype  = 'fir';
            cfgtemp.lpfreq      = config{irat}.LFP.lpfilter_wod_detection;
            LFP_trial_lpfilt      = ft_preprocessing(cfgtemp, LFP_trial);
           
            
            %recover trial real timings to use it with muse markers
            starttrial              = LFP_trial.trialinfo.begsample / LFP_trial.fsample;
            endtrial                = LFP_trial.trialinfo.endsample / LFP_trial.fsample;
            offsettrial             = LFP_trial.trialinfo.offset / LFP_trial.fsample;
            
                          
                       
            %get hand-annotated wod timing
            wod_marker = MuseStruct{1}.markers.WOD.synctime(itrial);
            wod_delay  = MuseStruct{1}.markers.WOD.synctime(itrial) - MuseStruct{1}.markers.Vent_Off.synctime(itrial);
            
            
            
            chan_name= string(config{irat}.LFP.origin_WoD{itrial})
            chan_idx= strcmp(chan_name,(LFP_trial.label));
            
            
            %select times where to search WOD peak
            t = LFP_trial_lpfilt.time{1};
            t_1 = t > (wod_marker + config{irat}.LFP.wod_toisearch(1) - starttrial + offsettrial);
            t_2 = t < (wod_marker + config{irat}.LFP.wod_toisearch(2) - starttrial + offsettrial);
            t_sel = t_1 & t_2;
            
            [v_peak_wod, t_peak_wod] = findpeaks(-LFP_trial_lpfilt.trial{1}(chan_idx,t_sel),t(t_sel),'NPeaks',1,'SortStr','descend','WidthReference','Halfheight');
            clear t t_1 t_2 t_sel

            LFP_trial.time{1}= LFP_trial.time{1} - t_peak_wod;
            
            
    
           
    


%% Plot overdraw by channel for all rats

if irat== 4
fig_ori = figure;
else
end

hold on
color= 'k';

p2 = plot(LFP_trial.time{1},LFP_trial.trial{1}(chan_idx,:),'k','LineWidth',1);
p2.Color(4)= 0.2;

axis tight

xlim([-5 5]);
xlabel('Time for WoD (s)', 'Fontsize',15);
ylabel('Voltage (µV)', 'Fontsize',15);

hold on


end %itrial    
end %irat


imagesavedir= '\\lexport\iss01.charpier\analyses\wod\images\overdraw'


if ~isfolder(imagesavedir)
    fprintf('Creating directory %s\n',imagesavedir);
    mkdir(imagesavedir)
end


set(fig_ori,'PaperOrientation','landscape');
set(fig_ori,'PaperUnits','normalized');
set(fig_ori,'PaperPosition', [0 0 1 1]);
print(fig_ori, '-dpdf', fullfile(imagesavedir,sprintf('AllRats_orichans.pdf')),'-r600');
print(fig_ori, '-dpng', fullfile(imagesavedir,sprintf('AllRats_orichans.png')),'-r600');
%savefig(fig,fullfile(imagesavedir,sprintf('AllRats_allchans')));
close all
