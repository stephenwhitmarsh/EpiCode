function [data] = preprocessing_eeg_emg(cfg,ipart,idir,force)
% load the data : only the channels set in setparams
% inverse the data
% notch
% rereferencing "average"
% AJOUTER LOAD DIFFERENT DES EMG
% AJOUTER LES PARAMETRES DANS SETPARAMS

datapath = fullfile(cfg.rawdir,[cfg.directorylist{ipart}{idir},'.TRC']);
fname = fullfile(cfg.datasavedir,[cfg.prefix,'data_preprocessed.mat']);

if exist(fname,'file') && force == false
    fprintf('*******************************\n');
    fprintf('** Loading data preprocessed **\n');
    fprintf('*******************************\n\n');
    load(fname,'data');
else
    
    if force == true
        fprintf('*************************************\n');
        fprintf('** Forced redoing of preprocessing **\n');
        fprintf('*************************************\n\n');
    else
        fprintf('*******************\n');
        fprintf('** Preprocessing **\n');
        fprintf('*******************\n\n');
    end
    
    
    % EEG preprocessing
    cfgtemp             = [];
    cfgtemp             = cfg.preproc_eeg;
    cfgtemp.dataset     = datapath;
    cfgtemp.bsfilter    = 'yes';
    cfgtemp.bsfreq      = [49 51];
    data_EEG            = ft_preprocessing(cfgtemp);
    
    
    % EMG preprocessing
    cfgtemp             = [];
    cfgtemp.channel     = cfg.preproc_emg;
    cfgtemp.dataset     = datapath;
    data_EMG            = ft_preprocessing(cfgtemp);
    
    
    %Appendata
    cfgtemp                 = [];
    cfgtemp.keepsampleinfo  = 'yes';
    data                    = ft_appenddata(cfgtemp, data_EEG, data_EMG);
    
    %inverse the data according to standard clinical visualisation
    data.trial{1} = cfg.inversedata*data.trial{1};
    
    
    % Plot to check if preprocessing was good :
    %         figure;
    %         hold;
    %
    %         offset1239 = 1239*256; %sample of beginning of the plot
    %         offset1269 = 1269*256; %sample of endning of the plot
    %         for i = 1:length(data.label)
    %         plot(data.time{1}(offset1239:offset1269),data.trial{1}(i,offset1239:offset1269)+(500*i));
    %         end
    %
    %         xlim([1239,1269]);
    %         ylabel('Channel name');
    %
    %         tick = 500;
    %         yticks(0 : tick : length(data.label)*500);
    %         set(gca, 'YTickLabel',[data.label(end)', data.label(1:end-1)']); %why Y axis not in the good order
    %         set(gca,'TickDir','out');
    
    
    % Save
    save(fname,'data');
    
end
end


