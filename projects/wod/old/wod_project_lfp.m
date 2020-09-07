function wod_project_lfp(slurm_task_id)

if ispc
    addpath (genpath('\\lexport\iss01.charpier\analyses\wod\EpiCode\shared'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\wod\EpiCode\external'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\wod\EpiCode\templates'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\wod\EpiCode\projects\wod'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\wod\EpiCode\projects\dtx'))
    addpath \\lexport\iss01.charpier\analyses\wod\fieldtrip-20200607
    addpath \\lexport\iss01.charpier\analyses\wod\CEDS64ML
    CEDS64LoadLib('\\lexport\iss01.charpier\analyses\wod\CEDS64ML');
    
elseif isunix
    addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/EpiCode/shared'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/EpiCode/external'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/EpiCode/templates'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/EpiCode/projects/wod'))
    addpath /network/lustre/iss01/charpier/analyses/wod/fieldtrip-20200607
end


ft_defaults

config = wod_setparams_lfp;


%% read all data of the whole file, for all rats
for irat = 1:size(config,2)
    
    %read all timing of events
    CEDStruct{irat}               = readCEDmarkers(config{irat}, false);

    for ipart = 1:size(config{irat}.directorylist,2) %list of protocols
        for idir = 1:size(config{irat}.directorylist{ipart},2) %nr of files per protocol. always 1 in this case
            data{irat}{ipart}{idir} = {};
            
            %find the data file. search for .smr or .smrx file
            temp        = dir(fullfile(config{irat}.rawdir, [config{irat}.directorylist{ipart}{idir}, '.smr*']));
            datapath    = fullfile(config{irat}.rawdir, temp.name);

            %open the file in Spike2.
            fid = CEDS64Open(datapath);
            if fid<0
                error('error while opening file %s. \nCheck that the path is correct, and that this file is not already opened in Spike2.',datapath);
            end
            
            for ichan = 1:CEDS64MaxChan(fid)
                chantitle = [];
                [~, chantitle] = CEDS64ChanTitle(fid, ichan);
                
                if strfind(chantitle, 'Power')
                    data{irat}{ipart}{idir}{end+1} = readCEDcontinous(config{irat}, chantitle, ipart, idir);
                end
                    
            end
        end
    end
end

%% save CEDStruct all rats and data all rats
 


%% load pre-computed data and analyse
%             
% 
% dat = readCEDcontinous(config{irat},config{irat}.LFP.channel{imarker}{ichannel}, ipart,idir);
% 
% for irat = 1
%     ipart = 1;
%     
%     %% Get LFP data
%     
%     % 	For readCEDMarkers and CED2mat : Can only be done on windows because
%     %   of CEDS64 library. For Linux, set force argument to false, it will
%     %   load the results.
%     
%     
%     %si besoin : ajouter ipart aux noms d'images et de data enregistrées
%     
%     %read all events
%     %directorylist CED
%     %datasavedir pour data coupées et event
%     
%     
%     
%     %demean
%     for ipart = 1:length(dat_LFP)
%         for imarker = 1:length(dat_LFP{ipart})
%             config{irat}temp                 = [];
%             config{irat}temp.demean          = config{irat}.LFP.demean{imarker};
%             config{irat}temp.baselinewindow  = config{irat}.LFP.baselinewindow{imarker};
%             dat_LFP{ipart}{imarker} = ft_preprocessing(config{irat}temp, dat_LFP{ipart}{imarker});
%         end
%     end
%     
%     %compute derivative of Vm
%     config{irat}temp             = [];
%     config{irat}temp.derivative  = 'yes';
%     der                 = ft_preprocessing(config{irat}temp, dat_LFP{1}{2});
%     
%     % plot data, ECoG-M1G
%     figure;
%     subplot(2,1,1);hold;
%     for itrial = 1:size(dat_LFP{1}{1}.trial,2)
%         plot(dat_LFP{1}{1}.time{itrial},-dat_LFP{1}{1}.trial{itrial});
%     end
%     xlim([-1 1]);
%     
%     
%     % plot data, enveloppe of Vm
%     subplot(2,1,2);hold;
%     for itrial = 1:size(dat_LFP{1}{1}.trial,2)
%         %        y = envelope(dat_LFP{1}{2}.trial{itrial},1000,'rms');
%         
%         y =  dat_LFP{1}{2}.trial{itrial};
%         t = dat_LFP{1}{2}.time{itrial};
%         [AP_vals, AP_locs] = findpeaks(y,t,'MinPeakProminence',20);
%         toremove = false(size(y));
%         for i = 1:size(AP_locs,2)
%             toremove = (toremove | (t>AP_locs(i)-0.001 & t<AP_locs(i)+0.001) );
%         end
%         y(toremove) = NaN;
%         toremove2 = der.trial{itrial}>0.2 | der.trial{itrial}<-0.2;
%         y(toremove2) = NaN;
%         
%         plot(dat_LFP{1}{2}.time{itrial},y);
%     end
%     xlim([-1 1]);
%     ylim([-10 40]);
%     
%     % sum intra
%     all_intra_summed = dat_LFP{1}{2}.trial{1};
%     for itrial = 2:size(dat_LFP{1}{2}.trial,2)
%         all_intra_summed = all_intra_summed + dat_LFP{1}{2}.trial{itrial};
%     end
%     plot(dat_LFP{1}{2}.time{1},all_intra_summed);
%     
%     %% Rm dyn
%     rmdyn = readtable('Z:\analyses\lgi1\DTX-INTRA\neurones-DTX\test2 RmDyn\Rm_dyn_vals.txt'); %output xy view from Spike2
%     rmdyn = table2array(rmdyn)';
%     %1st line : time
%     %2nd line : rm values
%     figure;
%     scatter(rmdyn(1,:), rmdyn(2,:),'.k');
%     rmdyn_movmean = movmean(rmdyn(2,:),20); %sur 200ms
%     plot(rmdyn(1,:), rmdyn_movmean,'LineWidth',2);
%     %     plot(rmdyn(1,:), rmdyn(2,:));
%     xlim([2416.004 2418.598])
%     
%     
    
end