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
%remove breathing and ekg channel
        cfgtemp         = [];
        cfgtemp.channel = {'all', '-E0', '-Respi', '-ECG','-Puff'};
        LFP             = ft_selectdata(cfgtemp, LFP);
        LFP_cleaned     = LFP; %save for later removing of artefacts        
        
        
        
end %irat

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
            
            LFP_trial.time{1}= LFP_trial.time{1} - wod_delay;
            
           
           
    




%% Plot overdraw by channel for all rats

% if irat == 4
% fig = figure;
% else
% end
% hold on
% color = 'k';%color of trials when plotted
% h = 10000;% FIXME à ajuster
% 
% 
% for ichan= 1 : size(LFP_trial.label,1)
%     chan_name = string(LFP_trial.label(ichan))
%     
%     chan_idx= strcmp(chan_name,(LFP_trial.label));
%     p = plot(LFP_trial.time{1}, LFP_trial.trial{1}(chan_idx, :)+(numel(LFP_trial.label)+1)*h-h*ichan,'k', 'LineWidth', 2);
% %     p.Color(4) = 0.2;
% end
% 
% 
% axis tight
% ylim([-h (numel(LFP_trial.label)+2)*h]);
% xlim([-10 20]);
% xlabel('Time for WoD (s)', 'Fontsize',15);
% ylabel('Channel name', 'Fontsize',15);
% tick = h;
% yticks(h : tick : numel(LFP_trial.label)*h);
% chan_name_all = string(LFP_trial.label)';
% set(gca, 'YTickLabel', flip(chan_name_all));
% set(gca, 'FontWeight','bold', 'Fontsize',15, 'TickDir','out', 'YTickLabel',flip(chan_name_all));
% 
% hold on

if irat== 4
fig_ori = figure;
else
end

hold on
color= 'k';

chan_name= string(config{irat}.LFP.origin_WoD{itrial})
chan_idx= strcmp(chan_name,(LFP_trial.label));
p2 = plot(LFP_trial.time{1},LFP_trial.trial{1}(chan_idx,:),'k','LineWidth',1);
p2.Color(4)= 0.2;

axis tight
ylim([-2000 500]);
xlim([-5 5]);
xlabel('Time for WoD (s)', 'Fontsize',15);
ylabel('Voltage (µV)', 'Fontsize',15);

hold on


end %itrial    
end %irat


imagesavedir= '\\lexport\iss01.charpier\analyses\wod\images\overdraw';


if ~isfolder(imagesavedir)
    fprintf('Creating directory %s\n',imagesavedir);
    mkdir(imagesavedir)
end

set(fig,'PaperOrientation','landscape');
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition', [0 0 1 1]);
print(fig, '-dpdf', fullfile(imagesavedir,sprintf('AllRats_allchans.pdf')),'-r600');
print(fig, '-dpng', fullfile(imagesavedir,sprintf('AllRats_allchans.png')),'-r600');

set(fig_ori,'PaperOrientation','landscape');
set(fig_ori,'PaperUnits','normalized');
set(fig_ori,'PaperPosition', [0 0 1 1]);
print(fig_ori, '-dpdf', fullfile(imagesavedir,sprintf('AllRats_orichans.pdf')),'-r600');
print(fig_ori, '-dpng', fullfile(imagesavedir,sprintf('AllRats_orichans.png')),'-r600');
savefig(fig,fullfile(imagesavedir,sprintf('AllRats_allchans')));
close all
