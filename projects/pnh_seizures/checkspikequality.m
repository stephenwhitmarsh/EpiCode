function [quality]  = checkspikequality(cfg,force)

warning('off','all')
ft_warning off

ldir = dir2(fullfile(cfg.eegdir,cfg.directory_searchstring));
ldir = ldir([ldir.isdir] == 1);
fname = fullfile(cfg.outputdir,[cfg.prefix,'_qualitycheck.mat']);

if exist(fname,'file') && force == false
    load(fname,'quality');
else
    
    quality = [];
    if exist(fname,'file')
        disp('Loading pre-exisiting data');
        load(fname,'quality');
        startindx = size(quality.fname,2);
        if startindx < 1
            startindx = 1;
        end
    else
        startindx = 1;
    end
    
    clear hasfiles
    for idir = 1 : size(ldir,1)
        
        % check if directory micro data
        if ~isempty(dir(fullfile(cfg.eegdir,ldir(idir).name,'*_m*.ncs')))
            hasfiles(idir) = 1;
        else
            fprintf('*** No micro data in %s\n',fullfile(cfg.eegdir,ldir(idir).name));
        end
    end
    
    ldir = ldir(find(hasfiles));
    
    % do not override previous calculations, or rather, only the last
    
    for idir = startindx : size(ldir,1)
        
        files = dir(fullfile(cfg.eegdir,ldir(idir).name,'*_m*.ncs'));
        
        try
            % load trials for selected MICRO channels
            for ifile = 1 : size(files,1)
                quality.fname{idir}{ifile}          = fullfile(cfg.eegdir,ldir(idir).name,files(ifile).name);
                fprintf('Loading: %s\n',quality.fname{idir}{ifile});
                
                quality.current.idir  = idir;
                quality.current.ifile = ifile;
                
                
                cfgtemp                             = [];
                cfgtemp.continuous                  = 'yes';
                cfgtemp.dataset                     = quality.fname{idir}{ifile};
                cfgtemp.hpfilter                    = 'yes';
                cfgtemp.hpfreq                      = 500;
                dat                                 = ft_preprocessing(cfgtemp);
                
                fprintf('Calculating RMS & MAD\n');
                quality.dat_rms{idir}{ifile}        = rms(dat.trial{1});
                quality.dat_mad{idir}{ifile}        = mad(dat.trial{1});
                
                for imad = 4 : 8 % start at mad of 4
                    fprintf('Thresholding MAD%d of %s\n',imad,quality.fname{idir}{ifile});
                    %                     whos('patdatlist');
                    %                     whos('dat');
                    indxl = abs(dat.trial{1}) >= quality.dat_mad{idir}{ifile}*imad;
                    
                    % only once for every threshold passing
                    indxl(1) = 0;
                    indx = find(indxl);
                    temp.indx{idir}{ifile}{imad-3} = find(abs(dat.trial{1}(indx-1)) < quality.dat_mad{idir}{ifile}*imad);
                end
                clear dat
                
            end % type
            
            for imad = 4 : 8 % start at mad of 4
                for ifile = 1 : size(files,1)
                    overlap = [];
                    for ifile2 = 1 : size(files,1)
                        if ifile ~= ifile2
                            [C,IA,IB] = intersect(temp.indx{idir}{ifile}{imad-3},temp.indx{idir}{ifile2}{imad-3});
                            overlap = [overlap; IA];
                        end
                    end
                    quality.dat_mad_cnt{idir}{ifile}(imad-3)        = size(temp.indx{idir}{ifile}{imad-3},2);
                    quality.dat_mad_cnt_clean{idir}{ifile}(imad-3)  = size(setdiff(temp.indx{idir}{ifile}{imad-3},overlap'),2);
                end
            end
            
            clear temp
            
            % save and continue with next directory.
            save(fname,'quality');
        catch
            sprintf('something went wrong with: %s',quality.fname{idir}{ifile});
        end
    end
end