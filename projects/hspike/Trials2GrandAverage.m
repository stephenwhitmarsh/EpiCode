function [dat, dat_hyp, trialinfo] = Trials2GrandAverage(cfg, force)

fname_out = fullfile(cfg{1}.datasavedir, 'Trials2GrandAverage.mat');

if exist(fname_out,'file') && force == false
    fprintf('Loading results spikeTrialStats_GrandAverage\n');
    
    % repeat to deal with load errors
    count = 0;
    err_count = 0;
    while count == err_count
        try
            load(fname_out);
        catch ME
            err_count = err_count + 1;
        end
        count = count + 1;
    end
    return
end

hyplabels   = ["REM",  "AWAKE", "PHASE_1", "PHASE_2", "PHASE_3"];

% organize all data in single struct array

iunit = 1;
clear dat
trialinfo = table;

for ipatient = 1 : size(cfg, 2)
    
    SpikeTrials_timelocked  = readSpikeTrials_MuseMarkers(cfg{ipatient});
    SpikeDensity_timelocked = spikeTrialDensity(cfg{ipatient});
    
    if isempty(SpikeDensity_timelocked)
        continue
    end
    
    for ipart = 1 : size(SpikeDensity_timelocked, 2)
        
        if isempty(SpikeDensity_timelocked{ipart})
            continue
        end
        
        for markername = string(fields(SpikeDensity_timelocked{ipart}.sdf_bar))'
            
            for ilabel = 1 : size(SpikeDensity_timelocked{ipart}.sdf_bar.(markername).label, 2)
                
                fprintf('Adding: patient %d, part %d, %s, %s\n', ipatient, ipart, markername, SpikeDensity_timelocked{ipart}.sdf_bar.(markername).label{ilabel});
                dat{iunit}                          = [];
                dat{iunit}.avg                      = SpikeDensity_timelocked{ipart}.sdf_lin.(markername).avg(ilabel, :);
                dat{iunit}.time                     = SpikeDensity_timelocked{ipart}.sdf_lin.(markername).time;
                dat{iunit}.label{1}                 = 'Cluster';
                dat{iunit}.dimord                   = 'chan_time';
                
                trialinfo.ipatient(iunit)           = ipatient;
                trialinfo.ipart(iunit)              = ipart;
                trialinfo.ilabel(iunit)             = ilabel;
                trialinfo.markername(iunit)         = markername;
                trialinfo.cluster_group{iunit}      = SpikeTrials_timelocked{ipart}.(markername).cluster_group{ilabel};
                if strfind(SpikeTrials_timelocked{ipart}.(markername).cluster_group{ilabel}, "good")
                    trialinfo.good(iunit)           = true;
                else
                    trialinfo.good(iunit)           = false;
                end
                
                % check if unit responds statistically 
                trialinfo.responsive(iunit)         = false;
                if isfield(SpikeDensity_timelocked{ipart}.stat.(markername){ilabel}, 'posclusters')
                    for ipos = 1 : size(SpikeDensity_timelocked{ipart}.stat.(markername){ilabel}.posclusters, 2)
                        if SpikeDensity_timelocked{ipart}.stat.(markername){ilabel}.posclusters(ipos).prob < 0.01
                            trialinfo.responsive(iunit) = true;
                        end
                    end
                end
                if isfield(SpikeDensity_timelocked{ipart}.stat.(markername){ilabel}, 'negclusters')
                    for ineg = 1 : size(SpikeDensity_timelocked{ipart}.stat.(markername){ilabel}.negclusters, 2)
                        if SpikeDensity_timelocked{ipart}.stat.(markername){ilabel}.negclusters(ineg).prob < 0.01
                            trialinfo.responsive(iunit) = true;
                        end
                    end
                end
                
                SpikeTrials_timelocked{ipart}.(markername).trialinfo.hyplabel(SpikeTrials_timelocked{ipart}.(markername).trialinfo.hyplabel == "NO_SCORE") = "AWAKE";
                
                for hyplabel = hyplabels
                    
                    trials = SpikeTrials_timelocked{ipart}.(markername).trialinfo.hyplabel == hyplabel;
                    dat_hyp.(hyplabel){iunit}.avg                      = squeeze(nanmean(SpikeDensity_timelocked{ipart}.sdf_lin.(markername).trial(trials, ilabel, :), 1))';
                    dat_hyp.(hyplabel){iunit}.time                     = SpikeDensity_timelocked{ipart}.sdf_lin.(markername).time;
                    dat_hyp.(hyplabel){iunit}.label{1}                 = 'Cluster';
                    dat_hyp.(hyplabel){iunit}.dimord                   = 'chan_time';
                    dat_hyp.(hyplabel){iunit}.trialinfo.ipatient       = ipatient;
                    dat_hyp.(hyplabel){iunit}.trialinfo.ipart          = ipart;
                    dat_hyp.(hyplabel){iunit}.trialinfo.markername     = markername;
                    dat_hyp.(hyplabel){iunit}.trialinfo.cluster_group  = SpikeTrials_timelocked{ipart}.(markername).cluster_group{ilabel};
                    
                end 
                iunit = iunit + 1;
            end
        end
    end
end

clear trialinfo_sum
trialinfo_sum = table;
i = 1;
for ipatient = unique(trialinfo.ipatient)'
    for ipart = 1 : 3
        trialinfo_sum.ipatient(i)           = ipatient;
        trialinfo_sum.ipart(i)              = ipart;
        trialinfo_sum.units(i)              = sum(trialinfo.ipatient == ipatient & trialinfo.ipart == ipart);
        trialinfo_sum.responsive_p(i)       = round((sum(trialinfo.ipatient == ipatient & trialinfo.ipart == ipart & trialinfo.responsive == true) / sum(trialinfo.ipatient == ipatient & trialinfo.ipart == ipart)) * 100, 0);  
        trialinfo_sum.sua(i)                = sum(trialinfo.ipatient == ipatient & trialinfo.ipart == ipart & trialinfo.good == true);
        trialinfo_sum.mua(i)                = sum(trialinfo.ipatient == ipatient & trialinfo.ipart == ipart & trialinfo.good == false);      
        trialinfo_sum.sua_responsive_p(i)   = round((sum(trialinfo.ipatient == ipatient & trialinfo.ipart == ipart & trialinfo.good == true & trialinfo.responsive == true) / sum(trialinfo.ipatient == ipatient & trialinfo.ipart == ipart & trialinfo.good == true)) * 100, 0);
        trialinfo_sum.mua_responsive_p(i)   = round((sum(trialinfo.ipatient == ipatient & trialinfo.ipart == ipart & trialinfo.good == false & trialinfo.responsive == true) / sum(trialinfo.ipatient == ipatient & trialinfo.ipart == ipart & trialinfo.good == false)) * 100, 0);  
        i = i + 1;
    end
end

fname = fullfile(cfg{1}.datasavedir, 'overview_GA_all');
writetable(trialinfo, fname);
fname = fullfile(cfg{1}.datasavedir, 'overview_GA_sum');
writetable(trialinfo_sum, fname);

save(fname_out, 'dat', 'dat_hyp', 'trialinfo', '-v7.3');
