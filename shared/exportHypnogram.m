function exportHypnogram(cfg)

% EXPORT_HYPNOGRAM searches through patient directly, finds hypnogram files and
% NeuraLynx directories that overlap with hypnogram files, aligns MicroMed
% with NeuraLynx data, and writes new Events.mrk files for Muse in all
% NeuraLynx directories that contain hypnogram data.
%
% Use as
%    exportHypnogram(cfg)
%
% Necessary fields:
%
% cfg.patientdir            = [path to patient data, one above /eeg and /eegmicromed]
% cfg.hyp.micromedchannel   = [micromed channel to use for alignment]
% cfg.imagesavedir          = [path to directory in which to save plots]
% cfg.prefix                = [string to prefix output figure, e.g. patient ID]
% cfg.hyp.backupdir         = [directory to backup Muse marker files]
%
% Example:
%
% cfg                       = [];
% cfg.patientdir            = '/network/lustre/iss01/epimicro/patients/raw/pat_02711_1193';
% cfg.hyp.micromedchannel   = 'F3p6';
% cfg.imagesavedir          = '/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/images/hspike';
% cfg.prefix                = 'P1-';
% cfg.hyp.backupdir         = '/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/markerbackup';
% export_hypnogram(cfg);
%
% Note:
% Names of markers that contain a space (' ') or minus ('-') will be
% replaced by an underscore ('_').
%
% Dependencies: writeMuseMarkers.m, dir2.m, recent FieldTrip version

% This file is part of EpiCode, see
% http://www.github.com/stephenwhitmarsh/EpiCode for documentation and details.
%
%    EpiCode is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    EpiCode is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with EpiCode. If not, see <http://www.gnu.org/licenses/>.

feature('DefaultCharacterSet', 'UTF8') %# for all Character support, or 'CP1252'

% list of all hypnogram files in MicroMed directory
micromed_hypnfilelist                   = dir2(fullfile(cfg.patientdir,'eegmicromed','*.hypn'));

% get some info about the 'real' time of recording the MicroMed data
for ifile = 1 : size(micromed_hypnfilelist,1)

    [~,name,~]                              = fileparts(micromed_hypnfilelist(ifile).name);
    fname                                   = fullfile(cfg.patientdir,'eegmicromed',[name,'.bni']);
    hdrMM                                   = ft_read_header(fullfile(cfg.patientdir,'eegmicromed',[name,'.TRC']));

    opts                                    = detectImportOptions(fname,'FileType','text','Delimiter','=');
    opts.VariableNames                      = {'key','value'};
    opts.Delimiter                          = {'='};
    opts.VariableTypes                      = {'char','char'};
    T2                                      = readtable(fname,opts);
    micromed_date{ifile}                    = datetime(cell2mat(T2.value(find(strcmp(T2.key,'Date')))),'Format','MM/dd/yyyy');
    micromed_time{ifile}                    = datetime(cell2mat(T2.value(find(strcmp(T2.key,'Time')))),'Format','HH/mm/ss');
    micromed_starttime{ifile}               = datetime([cell2mat(T2.value(find(strcmp(T2.key,'Date')))), cell2mat(T2.value(find(strcmp(T2.key,'Time'))))],'Format','MM/dd/yyyyHH:mm:ss');
    micromed_starttime{ifile}               = datetime(micromed_starttime{ifile},'Format','yyy/MM/dd HH:mm:ss');
    micromed_endtime{ifile}                 = micromed_starttime{ifile} + seconds(hdrMM.nSamples / hdrMM.Fs);
end

% overview and timing of neuralynx data
neuralynx_dirlist = dir2(fullfile(cfg.patientdir,'eeg'));
neuralynx_dirlist = neuralynx_dirlist([neuralynx_dirlist.isdir]);

% add "real time" of computer clock to neuralynx file time
for idir = 1 : size(neuralynx_dirlist,1)
    neuralynx_datafiles{idir}   = dir2(fullfile(neuralynx_dirlist(idir).folder,neuralynx_dirlist(idir).name,'*m*.ncs'));
    neuralynx_txtfiles{idir}    = dir2(fullfile(neuralynx_dirlist(idir).folder,neuralynx_dirlist(idir).name,'*.txt'));
    neuralynx_hdr{idir}         = ft_read_header(fullfile(neuralynx_dirlist(idir).folder,neuralynx_dirlist(idir).name,neuralynx_datafiles{idir}(1).name)); % read header from first ncs file

    f = fopen(fullfile(neuralynx_txtfiles{idir}(1).folder,neuralynx_txtfiles{idir}(1).name));
    clear timestring
    while 1
        tline = fgetl(f);
        if ~ischar(tline), break, end
        searchstring = '## Time Opened (m/d/y)';
        try
            if strcmp(tline(1:length(searchstring)),searchstring)
                neuralynx_timestring = tline;
                disp('Great, found timestamp in header file');
                break
            end
        catch
            disp('Warning: something weird happened reading the txt time');
        end
    end
    fclose(f);

    neuralynx_timestring         = strsplit(neuralynx_timestring);
    neuralynx_headerdate{idir}   = [cell2mat(neuralynx_timestring(5)) ' ' cell2mat(neuralynx_timestring(7))];
    neuralynx_starttime{idir}    = datetime(neuralynx_headerdate{idir},'Format','MM/dd/yy HH:mm:ss.SSS');
    neuralynx_endtime{idir}      = neuralynx_starttime{idir} + seconds(neuralynx_hdr{idir}.nSamples /neuralynx_hdr{idir}.Fs);
end

% find overlap between micromed and neuralynx files
for iMicroMed = 1 : size(micromed_hypnfilelist,1)
    hasoverlap{iMicroMed} = [];

    % Search all neuralynx files
    for iNeuraLynx = 1 : size(neuralynx_datafiles,2)
        % look for any marker in the MicroMed marker file that corresponds
        % with the approximate time recorded in the NeuraLynx datafile
        try
            if isbetween(micromed_starttime{iMicroMed},neuralynx_starttime{iNeuraLynx},neuralynx_endtime{iNeuraLynx}) || isbetween(micromed_endtime{iMicroMed},neuralynx_starttime{iNeuraLynx},neuralynx_endtime{iNeuraLynx})
                hasoverlap{iMicroMed} = [hasoverlap{iMicroMed}, iNeuraLynx];
            end
        catch
        end
    end

    % check if there is a neuralynx file for every hypnogram file
    if isempty(hasoverlap{iMicroMed})
        fprintf('OOPS! Could not find Neuralynx file that overlaps with %s',micromed_hypnfilelist(iMicroMed).name);
    end
end

% loop over all micromed files
for iMicroMed = 1 : size(micromed_hypnfilelist,1)

    [~,name,~]      = fileparts(micromed_hypnfilelist(iMicroMed).name);
    cfgtemp         = [];
    cfgtemp.dataset = fullfile(micromed_hypnfilelist(iMicroMed).folder,[name,'.TRC']);
    cfgtemp.channel = cfg.hyp.micromedchannel;
    dat_MM          = ft_preprocessing(cfgtemp);

    dat_NL_concat   = [];
    MM_title        = [cfgtemp.dataset, ': ', cfgtemp.channel];
    NL_title        = [];

    % loop over (multiple) overlapping neuralynx files
    for idir = [hasoverlap{iMicroMed}]

        % figure out the same channel name for neuralynx
        channelnr = 0;
        l = dir2(fullfile(neuralynx_datafiles{idir}(iMicroMed).folder,'*.ncs'));

        for i = 1 : size(l,1)
            if findstr(l(i).name,['_',dat_MM.label{1}(1:end-1), '_', dat_MM.label{1}(end)])
                channelnr = i;
                fprintf('Found match: %s : %s\n',l(i).name,['_',dat_MM.label{1}(1:end-1), '_', dat_MM.label{1}(end)]);
            end
        end
        if channelnr == 0
            fprintf('Could not find similar channel\n');
        end

        cfgtemp             = [];
        cfgtemp.dataset     = fullfile(neuralynx_datafiles{idir}(1).folder,l(channelnr).name);
        cfgtemp.channel     = 'all';
        dat_NL              = ft_preprocessing(cfgtemp);
        if ~isempty(NL_title)
            NL_title = [NL_title, '\n'];
        end
        NL_title            = [NL_title, l(channelnr).name];

        % resample neuralynx to same samplerate of micromed data
        cfgtemp             = [];
        cfgtemp.resamplefs  = hdrMM.Fs;
        dat_NL              = ft_resampledata(cfgtemp,dat_NL);
        dat_NL_concat       = [dat_NL_concat dat_NL.trial{1}];
    end

    % neuralynx data is inverted from micromed
    dat_NL_concat = -dat_NL_concat;
    D = finddelay(dat_MM.trial{1},dat_NL_concat);

    X1 = 1:size(dat_MM.trial{1},2);
    X2 = (1:size(dat_NL_concat,2)) - D;

    fig = figure;
    fig.Renderer = 'Painters';
    ax1 = subplot(2,1,1);
    plot(X1,dat_MM.trial{1});
    title(MM_title,'Interpreter','none');
    ax = gca;
    ax.FontSize = 8;
    axis tight

    ax2 = subplot(2,1,2);
    plot(X2,dat_NL_concat);
    title(sprintf(NL_title),'Interpreter','none');
    ax = gca;
    ax.FontSize = 8;
    axis tight

    linkaxes([ax1,ax2]);
    set(fig,'PaperOrientation','landscape');
    set(fig,'PaperUnits','normalized');
    set(fig,'PaperPosition', [0 0 1 1]);

    [~,name,~] = fileparts(micromed_hypnfilelist(iMicroMed).name);
    print(fig, '-dpdf', fullfile(cfg.hyp.imagesavedir,[cfg.prefix,'alignment_',name,'.pdf']));
    print(fig, '-dpng', fullfile(cfg.hyp.imagesavedir,[cfg.prefix,'alignment_',name,'.png']),'-r600');
    close all

    % read hypnogram exported from Micromed
    hyp = readtable(fullfile(micromed_hypnfilelist(iMicroMed).folder,micromed_hypnfilelist(iMicroMed).name),'FileType','text');

    % apply synchronization (overwrite startSec and endSec)
    hyp.startSec    = round(hyp.startSec + (D / hdrMM.Fs),4);
    hyp.endSec      = round(hyp.endSec   + (D / hdrMM.Fs),4);

    % take the Neurlynx data headerinfo from the first file to know how many samples the file is
    filelength = 0;
    for idir = [hasoverlap{iMicroMed}]
        flist = dir2(fullfile(neuralynx_dirlist(idir).folder,neuralynx_dirlist(idir).name,'*.ncs'));
        datafile = fullfile(flist(1).folder,flist(1).name);
        hdr = ft_read_header(datafile);
        filelength = [filelength, hdr.nSamples / hdr.Fs];
    end

    % figure out which hypnogram markers belong to which Neuralynx file
    filelengthc     = cumsum(filelength);
    hyp.startfilenr = zeros(height(hyp),1);
    hyp.endfilenr   = zeros(height(hyp),1);
    for i = 1 : size(filelength,2)
        hyp.startfilenr(hyp.startSec >= filelengthc(i)) = hyp.startfilenr(hyp.startSec >= filelengthc(i)) + 1;
        hyp.endfilenr(hyp.endSec >= filelengthc(i))     = hyp.endfilenr(hyp.endSec >= filelengthc(i)) + 1;
    end

    % split those hypnogram markers span over multiple files
    if any(hyp.startfilenr ~= hyp.endfilenr)
        for imarker = find(hyp.startfilenr ~= hyp.endfilenr)
            hyp(end+1,:)            = hyp(imarker,:);
            hyp.startSec(end)       = filelengthc(hyp.endfilenr(imarker));
            hyp.startfilenr(end)    = hyp.endfilenr(end);
            hyp.endSec(imarker)     = filelengthc(hyp.endfilenr(imarker));
            hyp.endfilenr(imarker)  = hyp.startfilenr(imarker);
        end
    end

    % adjust timing of those hypnogram periods falling in multiple files
    for i = 1 : height(hyp)
        hyp.startSec(i) = hyp.startSec(i) - filelengthc(hyp.startfilenr(i));
        hyp.endSec(i)   = hyp.endSec(i) - filelengthc(hyp.startfilenr(i));
    end

    hypi = 1;

    % loop over relevant neuralynx directories
    for idir = [hasoverlap{iMicroMed}]

        %  read Muse events file
        name_mrk    = fullfile(neuralynx_dirlist(idir).folder,neuralynx_dirlist(idir).name,'Events.mrk');
        MuseStruct  = readMuseMarker(name_mrk);

        % select those hypnogram markers that belong to the current Muse (neuralynx) marker file
        hyp_file    = hyp(hyp.startfilenr == hypi,:);

        for label = unique(hyp_file.stage)'
            disp(['Working on ',label{1}]);

            if strcmp(label,'BEGIN')
                MuseStruct.markers.StartHypnogram.synctime       = hyp_file.startSec(strcmp(hyp_file.stage, label))';     % replace space with underscore
                MuseStruct.markers.StartHypnogram.trialnum       = 0;
                MuseStruct.markers.StartHypnogram.classgroupid   = '+3';
                MuseStruct.markers.StartHypnogram.comment        = 'Exported from hypnogram';
                MuseStruct.markers.StartHypnogram.editable       = 'Yes';
                MuseStruct.markers.StartHypnogram.classid        = '+666'; % will be replaced by writeMuseMarkers.m
                MuseStruct.markers.StartHypnogram.color          = 'black';
            end

            if strcmp(label,'END')

                MuseStruct.markers.EndHypnogram.synctime        = hyp_file.startSec(strcmp(hyp_file.stage, label))';     % replace space with underscore
                MuseStruct.markers.EndHypnogram.trialnum        = 0;
                MuseStruct.markers.EndHypnogram.classgroupid    = '+3';
                MuseStruct.markers.EndHypnogram.comment         = 'Exported from hypnogram';
                MuseStruct.markers.EndHypnogram.editable        = 'Yes';
                MuseStruct.markers.EndHypnogram.classid         = '+666'; % will be replaced by writeMuseMarkers.m
                MuseStruct.markers.EndHypnogram.color           = 'black';
            end


            % for those markers that have a duration
            if hyp_file.startSec(strcmp(hyp_file.stage, label)) ~= hyp_file.endSec(strcmp(hyp_file.stage, label))

                % replace markers in markerfile
                MuseStruct.markers.([strrep(label{1},' ','_'), '__START__']).synctime = [];
                MuseStruct.markers.([strrep(label{1},' ','_'), '__END__']).synctime = [];

                try
                    MuseStruct.markers.([strrep(label{1},' ','_'), '__START__']).synctime       = [MuseStruct.markers.([strrep(label{1},' ','_'), '__START__']).synctime, hyp_file.startSec(strcmp(hyp_file.stage, label))'];     % replace space with underscore
                catch
                    MuseStruct.markers.([strrep(label{1},' ','_'), '__START__']).synctime       = hyp_file.startSec(strcmp(hyp_file.stage, label))';     % replace space with underscore
                end
                MuseStruct.markers.([strrep(label{1},' ','_'), '__START__']).classgroupid   = '+3';
                MuseStruct.markers.([strrep(label{1},' ','_'), '__START__']).comment        = 'Exported from hypnogram';
                MuseStruct.markers.([strrep(label{1},' ','_'), '__START__']).editable       = 'Yes';
                MuseStruct.markers.([strrep(label{1},' ','_'), '__START__']).classid        = '+666'; % will be replaced by writeMuseMarkers.m

                switch label{1}
                    case 'AWAKE'
                        MuseStruct.markers.([strrep(label{1},' ','_'), '__START__']).color = 'cyan';
                    case 'PHASE 1'
                        MuseStruct.markers.([strrep(label{1},' ','_'), '__START__']).color = 'green';
                    case 'PHASE 2'
                        MuseStruct.markers.([strrep(label{1},' ','_'), '__START__']).color = 'blue';
                    case 'PHASE 3'
                        MuseStruct.markers.([strrep(label{1},' ','_'), '__START__']).color = 'red';
                    case 'PHASE 4'
                        MuseStruct.markers.([strrep(label{1},' ','_'), '__START__']).color = 'purple';
                    otherwise
                        MuseStruct.markers.([strrep(label{1},' ','_'), '__START__']).color = 'black';
                end

                % concatinates with existing markers. Duplicates are
                % removed by writeMuseMarker
                try
                    MuseStruct.markers.([strrep(label{1},' ','_'), '__END__']).synctime     = [MuseStruct.markers.([strrep(label{1},' ','_'), '__END__']).synctime, hyp_file.endSec(strcmp(hyp_file.stage, label))'];     % replace space with underscore
                catch
                    MuseStruct.markers.([strrep(label{1},' ','_'), '__END__']).synctime     = hyp_file.endSec(strcmp(hyp_file.stage, label))';     % replace space with underscore
                end
                MuseStruct.markers.([strrep(label{1},' ','_'), '__END__']).classgroupid = '+3';
                MuseStruct.markers.([strrep(label{1},' ','_'), '__END__']).comment      = 'Exported from hypnogram';
                MuseStruct.markers.([strrep(label{1},' ','_'), '__END__']).editable     = 'Yes';
                MuseStruct.markers.([strrep(label{1},' ','_'), '__END__']).classid      = '+666'; % will be replaced by writeMuseMarkers.m

                switch label{1}
                    case 'AWAKE'
                        MuseStruct.markers.([strrep(label{1},' ','_'), '__END__']).color = 'yellow';
                    case 'PHASE 1'
                        MuseStruct.markers.([strrep(label{1},' ','_'), '__END__']).color = 'green';
                    case 'PHASE 2'
                        MuseStruct.markers.([strrep(label{1},' ','_'), '__END__']).color = 'blue';
                    case 'PHASE 3'
                        MuseStruct.markers.([strrep(label{1},' ','_'), '__END__']).color = 'red';
                    case 'PHASE 4'
                        MuseStruct.markers.([strrep(label{1},' ','_'), '__END__']).color = 'purple';
                    otherwise
                        MuseStruct.markers.([strrep(label{1},' ','_'), '__END__']).color = 'black';
                end

                % round to 0 decimals
                MuseStruct.markers.([strrep(label{1},' ','_'), '__START__']).synctime  =  round(MuseStruct.markers.([strrep(label{1},' ','_'), '__START__']).synctime);
                MuseStruct.markers.([strrep(label{1},' ','_'), '__END__']).synctime    =  round(MuseStruct.markers.([strrep(label{1},' ','_'), '__END__']).synctime);

                % remove duplicates
                [~,IA,IC] = unique([MuseStruct.markers.([strrep(label{1},' ','_'), '__START__']).synctime; MuseStruct.markers.([strrep(label{1},' ','_'), '__END__']).synctime]','rows');
                MuseStruct.markers.([strrep(label{1},' ','_'), '__START__']).synctime   = MuseStruct.markers.([strrep(label{1},' ','_'), '__START__']).synctime(IA);
                MuseStruct.markers.([strrep(label{1},' ','_'), '__END__']).synctime     = MuseStruct.markers.([strrep(label{1},' ','_'), '__END__']).synctime(IA);

                % add trialnr (zeros)
                MuseStruct.markers.([strrep(label{1},' ','_'), '__START__']).trialnum   = zeros(size(MuseStruct.markers.([strrep(label{1},' ','_'), '__START__']).synctime))';
                MuseStruct.markers.([strrep(label{1},' ','_'), '__END__']).trialnum     = zeros(size(MuseStruct.markers.([strrep(label{1},' ','_'), '__END__']).synctime))';  % replace space with underscore

%             else
%                 MuseStruct.markers.(strrep(label{1},' ','_')).synctime      = hyp_file.endSec(strcmp(hyp_file.stage, label));  % replace space with underscore
%                 MuseStruct.markers.(strrep(label{1},' ','_')).trialnum      = zeros(size(MuseStruct.markers.(strrep(label{1},' ','_')).synctime));  % replace space with underscore
%                 MuseStruct.markers.(strrep(label{1},' ','_')).classgroupid  = '+3';
%                 MuseStruct.markers.(strrep(label{1},' ','_')).comment       = 'Exported from hypnogram';
%                 MuseStruct.markers.(strrep(label{1},' ','_')).editable      = 'Yes';
%                 MuseStruct.markers.(strrep(label{1},' ','_')).classid       = '+666'; % will be replaced by writeMuseMarkers.m
%                 switch label{1}
%                     case 'AWAKE'
%                         MuseStruct.markers.(strrep(label{1},' ','_')).color = 'yellow';
%                     case 'PHASE 1'
%                         MuseStruct.markers.(strrep(label{1},' ','_')).color = 'green';
%                     case 'PHASE 2'
%                         MuseStruct.markers.(strrep(label{1},' ','_')).color = 'blue';
%                     case 'PHASE 3'
%                         MuseStruct.markers.(strrep(label{1},' ','_')).color = 'red';
%                     case 'PHASE 4'
%                         MuseStruct.markers.(strrep(label{1},' ','_')).color = 'purple';
%                     otherwise
%                         MuseStruct.markers.(strrep(label{1},' ','_')).color = 'black';
%                 end



            end
        end

        % backup markerfile
        if ~exist(cfg.hyp.backupdir,'DIR')
            error('Backup directory does not exist');
        end
        if ~exist(fullfile(cfg.hyp.backupdir,neuralynx_dirlist(idir).name),'DIR')
            fprintf('Creating directory: %s\n',fullfile(cfg.hyp.backupdir,neuralynx_dirlist(idir).name));
            eval(sprintf('!mkdir %s',fullfile(cfg.hyp.backupdir,neuralynx_dirlist(idir).name)));
        end
        fname_backup = sprintf('Events_%s.mrk', datestr(now, 'mm-dd-yyyy_HH-MM-SS'));
        eval(sprintf('!cp %s %s',name_mrk,fullfile(cfg.hyp.backupdir,neuralynx_dirlist(idir).name,fname_backup)));
        fprintf('Succesfully backed up markerfile to %s\n',fullfile(cfg.hyp.backupdir,neuralynx_dirlist(idir).name,fname_backup));

        % write new event structure
        fname_out = fullfile(neuralynx_dirlist(idir).folder,neuralynx_dirlist(idir).name,'Events.mrk');
%         fname_out = fullfile(cfg.hyp.markerdir,[neuralynx_dirlist(idir).name,'.mrk']);

        writeMuseMarkerfile(MuseStruct, fname_out);
        hypi = hypi + 1;
    end

end % micromed files
