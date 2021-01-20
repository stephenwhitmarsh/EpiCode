
warning('off','all')
addpath /network/lustre/iss01/charpier/stephen.whitmarsh/WhitmarshEpilepsy
datapath = '/network/lustre/iss01/epimicro/patients/raw';
outputdir = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/qualitycheck';

% loop through patients
clear patlist patdatlist prefix
i = 1;
dlist = dir2(fullfile(datapath,'pat*'));
for idlist = 1 : size(dlist,1)
    fprintf('Reading %d of %d\n',idlist,size(dlist,1));
    % if EEG directory exist
    if isdir(fullfile(dlist(idlist).folder,dlist(idlist).name,'eeg'))
        
        slist = dir2(fullfile(dlist(idlist).folder,dlist(idlist).name,'eeg'));
        ii = 1;
        hasmicro = false;
        
        % loop through data directories
        for islist = 1 : size(slist,1)
            
            % if it is a directory dive in
            if isdir(fullfile(slist(islist).folder,slist(islist).name))
                
                % look for micro files
                flist = dir2(fullfile(slist(islist).folder,slist(islist).name,'*_m*.ncs'));
                if ~isempty(flist)
                    patdatlist{i}{ii} = slist(islist);
                    ii = ii + 1;
                    
                    
                    % found one so add it to the list
                    hasmicro = true;
                end
            end
        end
        
        if hasmicro
            % extract name for prefix
            prefix{i} =  dlist(idlist).name(1:9);
            
            % go to next entry
            i = i + 1;
        end
        
    end
end

% ssh stephen.whitmarsh@login01 ?
% alias sacct="sacct --units=G -o submit,jobid,partition,NodeList,user,jobname%25,reqmem,timelimit,elapsed,maxrss,averss,state,start"
% salloc -p normal -n 1 -c 28 --mem=120G -t 999:00:00 -J "qualitycheck"
% squeue -u $USER
% ssh lmb021
% cd /network/lustre/iss01/charpier/stephen.whitmarsh/WhitmarshEpilepsy/
% module load MATLAB
% matlab < plotspikequality.m 



maxNumCompThreads(1); % should be one, probably
parpool(28) % should be same as nr. of cores per node. 

parfor i = 1 : size(prefix,2)
    display(['Running parfor loop ', num2str(i)])
    cfg                         = [];
    cfg.eegdir                  = patdatlist{i}{1}.folder;
    cfg.outputdir               = outputdir;
    cfg.prefix                  = prefix{i};
    cfg.directory_searchstring  = '*';
    quality{i}                  = checkspikequality(cfg,1);
    
end

% 2614 = 37


% 
% 
% cfg                         = [];
% cfg.eegdir                  = '/network/lustre/iss01/epimicro/patients/raw/pat_02476_0929/eeg';
% cfg.outputdir               = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/';
% cfg.prefix                  = 'pat2476-';
% cfg.directory_searchstring  = '*';
% quality                     = checkspikequality(cfg,1);
% 
% 
% i = 1;
% for idir = 1 : size(quality.dat_mad_cnt,2)
%     for imad = 1 : 10
%         for ifile = 1 : size(quality.dat_mad_cnt{idir},2)
%             q(imad,ifile,idir) = quality.dat_mad_cnt{idir}{ifile}(imad);
% 
%         end
%     end
% end
% 
% q_norm = q ./ max(q,[],1);
% q_norm = q ;
% 
% labels = quality.fname{1};
% for i = 1 : size(labels,2)
%     indx = findstr(labels{i},'_m');
%     labels{i} = labels{i}(indx+1:end-4);
%     labels{i} = strrep(labels{i},'_','');
% end
% 
% for i = 1 : size(quality.fname,2)
%     temp = regexp(quality.fname{i}{1},filesep,'split');
%     temp = strrep(temp,'_',' ');
%     fnames{i} = temp{end-1}
% end
% 
% 
% fig = figure; 
% 
% qsel = squeeze(q(6,:,:));
% 
% subplot(3,1,1);
% imagesc(log(qsel))
% yticks(1:size(qsel,1))
% yticklabels(labels);
% xticklabels([]);
% 
% axis tight
% set(gca,'fontsize',6);
% title('Max-Min all');
% 
% subplot(3,1,2);
% imagesc(qsel./max(qsel,[],1))
% yticks(1:size(qsel,1))
% yticklabels(labels);
% xticklabels([]);
% 
% axis tight
% set(gca,'fontsize',6);
% title('Scaled over time');
% 
% subplot(3,1,3);
% imagesc(qsel./max(qsel,[],2))
% yticks(1:size(qsel,1))
% yticklabels(labels);
% xticks(1:size(qsel,2))
% xticklabels(fnames);
% xtickangle(90)
% axis tight
% set(gca,'fontsize',6);
% title('Scaled over channels');
% 
% fname_out = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/qualitycheck.pdf';
% 
% set(fig,'PaperOrientation','landscape');
% set(fig,'PaperUnits','normalized');
% set(fig,'PaperPosition', [0 0 1 1]);
% print(fig, '-dpdf', fullfile(fname_out));