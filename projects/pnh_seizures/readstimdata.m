function [markerdat_micro, markerdat_macro] = readstimdata(cfg,micromed_markers)

fname_output = fullfile(cfg.datasavedir,['stimdat_',cfg.postfix,'.mat']);
if exist(fname_output,'file') && cfg.force == false
    fprintf('Loading precomputed data: %s\n',fname_output);
    load(fname_output,'markerdat_micro','markerdat_macro');
else
    
    % loop over stimulation events
    for imarker = 1 : height(micromed_markers)
        
        flist_NL = dir2(fullfile(micromed_markers.DataFileFolder_NL(imarker,:),'*.ncs'));
        
        icounter = 1;
        for ifile = 1 : size(flist_NL,1)
            if ~contains(flist_NL(ifile).name,'SYNC') && ~contains(flist_NL(ifile).name,'SLI')
                
                % single electrode from directory (time) corresponding to marker
                fname           = fullfile(micromed_markers.DataFileFolder_NL(imarker,:),flist_NL(ifile).name);
                hdr             = ft_read_header(fname);
                trl_prestim     = cfg.prestim * hdr.Fs;
                trl_poststim    = cfg.poststim * hdr.Fs;
                
                if (micromed_markers.(cfg.samplecolumnNL)(imarker) + trl_poststim) < hdr.nSamples && (micromed_markers.(cfg.samplecolumnNL)(imarker) - trl_prestim) > 0
                    cfgtemp             = [];
                    cfgtemp.dataset     = fname;
                    cfgtemp.trl         = [micromed_markers.(cfg.samplecolumnNL)(imarker) - trl_prestim, micromed_markers.(cfg.samplecolumnNL)(imarker) + trl_poststim, -trl_prestim];
                    
                    if cfgtemp.trl(1,2) < 1
                        cfgtemp.trl(1) = 1;
                    end
                    if cfgtemp.trl(1,2) > hdr.nSamples
                        cfgtemp.trl(2) = hdr.nSamples;
                    end
                    
                    fprintf('Loading file %d of %d, for marker %d of %d',ifile, size(flist_NL,1), imarker, height(micromed_markers));
                    dat_file{icounter}  = ft_preprocessing(cfgtemp);
                else
                    fprintf('Stimulation extends beyond file!');
                end
                icounter = icounter + 1;
            end
        end
        
        % combine all channels for single stimulation event
        markerdat_micro{imarker} = ft_appenddata([],dat_file{:});
        markerdat_micro{imarker}.fsample = hdr.Fs;
        clear dat_file
        
        % read macro
        fname           = fullfile(micromed_markers.DateFileFolder(imarker,:),[micromed_markers.DataFileName(imarker,1:end-3),'TRC']);
        hdr             = ft_read_header(fname);
        trl_prestim     = cfg.prestim   * hdr.Fs;
        trl_poststim    = cfg.poststim  * hdr.Fs;
        
        cfgtemp         = [];
        cfgtemp.dataset = fname;
        cfgtemp.trl     = [micromed_markers.(cfg.samplecolumnMM)(imarker) - trl_prestim, micromed_markers.(cfg.samplecolumnMM)(imarker) + trl_poststim, -trl_prestim];
          
        if cfgtemp.trl(1,2) < 1
            cfgtemp.trl(1) = 1;
        end
        if cfgtemp.trl(1,2) > hdr.nSamples
            cfgtemp.trl(2) = hdr.nSamples;
        end
        
        markerdat_macro{imarker} = ft_preprocessing(cfgtemp);
        markerdat_macro{imarker}.fsample = hdr.Fs;
        
    end
    
    % remove date from label to concatinate
    for imarker = 1 : height(micromed_markers)
        for ilabel = 1 : size(markerdat_micro{imarker}.label,1)
            markerdat_micro{imarker}.label{ilabel} = markerdat_micro{imarker}.label{ilabel}(end-6:end);
        end
    end
    
    save(fname_output,'markerdat_micro','markerdat_macro','-v7.3');
    
end
fprintf('Done.\n');
