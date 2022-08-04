function Figure_PSTH_all


restoredefaultpath
if isunix
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/git/fieldtrip
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/projects/hspike/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/shared/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/shared/utilities
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/subaxis
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/fieldtrip
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/sigstar-master
    addpath(genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/scripts/releaseDec2015/'));
    addpath(genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/scripts/epishare-master'));
    addpath(genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/scripts/SPIKY_apr_2021'))
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/cbrewer/cbrewer
    
end

if ispc
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\git\fieldtrip
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\projects\hspike
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\shared
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\shared\utilities
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external\altmany-export_fig-8b0ba13\
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external\subaxis
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external\fieldtrip
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external\sigstar-master
    addpath(genpath('\\lexport\iss01.charpier\analyses\stephen.whitmarsh\scripts\epishare-master'));
    addpath(genpath('\\lexport\iss01.charpier\analyses\stephen.whitmarsh\scripts\SPIKY_apr_2021'));
    addpath          \\lexport\iss01.charpier\analyses\stephen.whitmarsh\scripts\MatlabImportExport_v6.0.0
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external\cbrewer\cbrewer
    
end

ft_defaults

config = hspike_setparams;

%% load data
for ipatient = 1 : 8
    config{ipatient}.spike.name  = ["template1", "template2", "template3", "template4", "template5", "template6"];
%     SpikeTrials{ipatient}        = readSpikeTrials(config{ipatient});
%     config{ipatient}.spike.name  = ["window"];
%     SpikeStats{ipatient}         = spikeTrialStats(config{ipatient});
    SpikePSTH{ipatient}       = spikePSTH(config{ipatient});
%     config{ipatient}.LFP.name    = ["template1", "template2", "template3", "template4", "template5", "template6"];
%     config{ipatient}.LFP.postfix = {'_all'};
%     LFPavg{ipatient}             = readLFPavg(config{ipatient});
%     SpikeWaveforms{ipatient}     = readSpikeWaveforms(config{ipatient});
end


%% average over templates and express as percentage change

for ipatient = 1 : 8
    for ipart = 1 : 3
        fn = fields(SpikePSTH{ipatient}{ipart}.psth);
        temp = [];
        for markername = string(fn)'
            if strcmp(markername, "all")
                continue
            end
            temp = [temp; permute(SpikePSTH{ipatient}{ipart}.psth.(markername).trial, [1,3,2])];
            t = SpikePSTH{ipatient}{ipart}.psth.(markername).time >= config{ipatient}.stats.bl(1).(markername)(1) & SpikePSTH{ipatient}{ipart}.psth.(markername).time <= config{ipatient}.stats.bl(1).(markername)(2);    
        end
        SpikePSTH{ipatient}{ipart}.psth.all.avg_norm = permute(nanmean(temp, 1), [3,2,1]);
        
        % select baseline
        bl = nanmean(SpikePSTH{ipatient}{ipart}.psth.all.avg_norm(:,t), 2);
        for iunit = 1 : size(bl, 1)
            SpikePSTH{ipatient}{ipart}.psth.all.avg_norm(iunit, :) = SpikePSTH{ipatient}{ipart}.psth.all.avg_norm(iunit, :) ./ bl(iunit);
%             SpikePSTH{ipatient}{ipart}.psth.all.avg_norm(iunit, :) = SpikePSTH{ipatient}{ipart}.psth.all.avg_norm(iunit, :) ./ max(abs(SpikePSTH{ipatient}{ipart}.psth.all.avg_norm(iunit, :)));
        end
    end
end

alpha = 0.001;

for ipatient = 1 : 8
    for ipart = 1 : 3
        
        fn = fields(SpikePSTH{ipatient}{ipart}.psth);       

        % determine responsiveness
        for iunit = 1 : size(SpikePSTH{ipatient}{ipart}.stat.(string(fn{1})), 2)

            SpikePSTH{ipatient}{ipart}.psth.all.responsive(iunit) = false;
            
            % if responsive in any template
            for markername = string(fn)'
                if strcmp(markername, "all")
                    continue
                end
                
                if isfield(SpikePSTH{ipatient}{ipart}.stat.(markername){iunit}, 'posclusters')
                    for ipos = 1 : size(SpikePSTH{ipatient}{ipart}.stat.(markername){iunit}.posclusters, 2)
                        if SpikePSTH{ipatient}{ipart}.stat.(markername){iunit}.posclusters(ipos).prob < alpha
                            SpikePSTH{ipatient}{ipart}.psth.all.responsive(iunit) = true;
                        end
                    end
                end
                if isfield(SpikePSTH{ipatient}{ipart}.stat.(markername){iunit}, 'negclusters')
                    for ineg = 1 : size(SpikePSTH{ipatient}{ipart}.stat.(markername){iunit}.negclusters, 2)
                        SpikePSTH{ipatient}{ipart}.stat.(markername){iunit}.negclusters(ineg).prob
                        if SpikePSTH{ipatient}{ipart}.stat.(markername){iunit}.negclusters(ineg).prob < alpha
                            SpikePSTH{ipatient}{ipart}.psth.all.responsive(iunit) = true;
                        end
                    end
                end   
            end
            
        end
    end
end
           
clear SpikePSTH_all
SpikePSTH_all.dat = [];
SpikePSTH_all.ipatient = [];
SpikePSTH_all.ipart = [];
SpikePSTH_all.iunit = [];
SpikePSTH_all.responsive = [];

for ipatient = 1 : 8
%     if ipatient == 7
%         continue
%     end
    
    for ipart = 1 : 3
        SpikePSTH_all.dat           = [SpikePSTH_all.dat;           SpikePSTH{ipatient}{ipart}.psth.all.avg_norm];
        SpikePSTH_all.ipatient      = [SpikePSTH_all.ipatient;      repmat(ipatient, [size(SpikePSTH{ipatient}{ipart}.psth.all.avg_norm, 1),1])];
        SpikePSTH_all.ipart         = [SpikePSTH_all.ipart;         repmat(ipart,    [size(SpikePSTH{ipatient}{ipart}.psth.all.avg_norm, 1),1])];
        SpikePSTH_all.responsive    = [SpikePSTH_all.responsive;    SpikePSTH{ipatient}{ipart}.psth.all.responsive'];       
    end
end

figure; 
subplot(1,4,1); image(SpikePSTH_all.dat, 'CDataMapping','scaled')
caxis([0, 4]);

subplot(1,4,2); imagesc(SpikePSTH_all.responsive)
subplot(1,4,3); imagesc(SpikePSTH_all.ipatient)
subplot(1,4,4); imagesc(SpikePSTH_all.ipart)

%%
figure
[X,Y] = meshgrid(-10:10);
Z = X + Y;
s = surf(X,Y,Z);
xlabel('X');
ylabel('Y');
zlabel('Z = C');
colorbar
 caxis([0 20])












SpikeMask_all = [];
SpikeMasknr_all = [];
alpha = 0.05;
i = 0;
for ipatient = 1 : 8
    for ipart = 1 : 3
        for markername = ["template1", "template2", "template3", "template4", "template5", "template6"]
            
            if ~isfield(SpikePSTH{ipatient}{ipart}.stat, markername)
                continue
            end
            if isempty(SpikePSTH{ipatient}{ipart}.stat.(markername))
                continue
            end
            i = i + 1;
            SpikeMasknr_all = [SpikeMasknr_all; i];

            for iunit = 1 : size(SpikePSTH{ipatient}{ipart}.stat.(markername), 2)
                
                mask = zeros(1, 121);
                
                if isfield(SpikePSTH{ipatient}{ipart}.stat.(markername){iunit}, 'posclusters')
                    for icluster = 1 : size(SpikePSTH{ipatient}{ipart}.stat.(markername){iunit}.posclusters, 2)
                        if SpikePSTH{ipatient}{ipart}.stat.(markername){iunit}.posclusters(icluster).prob < alpha
                            mask(SpikePSTH{ipatient}{ipart}.stat.(markername){iunit}.posclusterslabelmat == icluster) = 1;
                        end
                    end
                end
                
                if isfield(SpikePSTH{ipatient}{ipart}.stat.(markername){iunit}, 'negclusters')
                    for icluster = 1 : size(SpikePSTH{ipatient}{ipart}.stat.(markername){iunit}.negclusters, 2)
                        if SpikePSTH{ipatient}{ipart}.stat.(markername){iunit}.negclusters(icluster).prob < alpha
                            mask(SpikePSTH{ipatient}{ipart}.stat.(markername){iunit}.negclusterslabelmat == icluster) = -1;
                        end
                    end
                end
                SpikeMask_all(i, :) = mask;
            end
        end
    end
end
figure; imagesc(SpikeMask_all);
