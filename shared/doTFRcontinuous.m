function [TFR] = doTFRcontinuous(cfg,MuseStruct,force)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [TFR] = doTFR(cfg,MuseStruct,force)
%
% (c) Stephen Whitmarsh (stephen.whitmarsh@gmail.com)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname_out = fullfile(cfg.datasavedir,[cfg.prefix,'TFR.mat']);

if exist(fname_out,'file') && force == false
    fprintf('************************************\n');
    fprintf('*** Loading precomputed TFR data ***\n');
    fprintf('************************************\n\n');
    load(fname_out,'TFR');
    
else
    fprintf('********************************\n');
    fprintf('*** (re-) computing TFR data ***\n');
    fprintf('********************************\n\n');
    
    for ipart = 1 : size(MuseStruct,2)
        
        % just define MuseStruct_macro here for only one part to simplify code
        % and similarity with no-parts
        fprintf('\n*** Working on part %d ***\n',ipart)
        
        % process channels separately
        chan_counter = 1;
        clear chandat
        
        for ichan = 1 : size(cfg.TFR.channel,2)
            
            clear dirdat
            
            % loop over all directories (time), concatinating channel
            for idir = 1 : size(MuseStruct{ipart},2)
                
                % find corresponding file
                d                         = dir2(fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir},['*',cfg.TFR.channel{ichan},'*.ncs']));
                
                % load data
                cfgtemp                   = [];
                cfgtemp.dataset           = fullfile(d.folder,d.name);
                fprintf('LOADING: %s\n',cfgtemp.dataset);
                clear fname
                fname{1}                  = cfgtemp.dataset;
                dirdat{idir}              = ft_read_neuralynx_interp(fname);
                
                % downsample data
                cfgtemp                   = [];
                cfgtemp.resamplefs        = 100;
                dirdat{idir}              = ft_resampledata(cfgtemp,dirdat{idir});
                
                % truncate label to make them equal over files
                dirdat{idir}.label{1}     = dirdat{idir}.label{1}(end-6:end); % can be replaced by circus.channel
                
                % save sampleinfo
                hdr  = ft_read_header(fullfile(d.folder,d.name));
                cfg.sampleinfo_TFR{ipart}{ichan}(idir,:) = [1 hdr.nSamples];
                
            end % idir
            
            % concatinate data over files
            chandat{chan_counter} = dirdat{1};
            for idir = 2 : length(MuseStruct{ipart})
                fprintf('Concatinating directory %d, channel %d\n',idir, ichan);
                chandat{chan_counter}.trial{1} = [chandat{chan_counter}.trial{1} dirdat{idir}.trial{1}];
                chandat{chan_counter}.time{1}  = [chandat{chan_counter}.time{1} (dirdat{idir}.time{1} + chandat{chan_counter}.time{1}(end))];
            end
            
            chan_counter = chan_counter + 1;
            
        end % ichan
        
        dat{ipart} = ft_appenddata([],chandat{:});
        
        % TFR
        % time frequency analysis
        cfgtemp                         = [];
        cfgtemp.channel                 = 'all'; %ichannel;
        cfgtemp.method                  = 'mtmconvol';
        cfgtemp.output                  = 'pow';
        cfgtemp.taper                   = 'hanning';
        cfgtemp.pad                     = 'nextpow2';
        cfgtemp.keeptrials              = 'yes';
        cfgtemp.foi                     = 1:1:30;
        cfgtemp.t_ftimwin               = ones(size(cfgtemp.foi))*60;
        cfgtemp.toi                     = 0:30:dat{ipart}.time{1}(end);
        TFR{ipart}                      = ft_freqanalysis(cfgtemp,dat{ipart});    
        
    end % ipart
    save(fname_out,'TFR','-v7.3');
    
end % if file already exists

