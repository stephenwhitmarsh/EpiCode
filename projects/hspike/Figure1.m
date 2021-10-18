function Figure1

MuseStruct{8}           = [];
LFP_cluster{8}          = [];
LFP_cluster_detected{8} = [];
clusterindx{8}          = [];

for ipatient = 1:8
    config                                                      = hspike_setparams;
    [clusterindx{ipatient}, LFP_cluster{ipatient}]              = clusterLFP(config{ipatient});
    [config{ipatient}, LFP_cluster{ipatient}]                   = alignClusters(config{ipatient},  LFP_cluster{ipatient}{1}.Hspike.kmedoids{6});
    [MuseStruct{ipatient}, ~, LFP_cluster_detected{ipatient}]   = detectTemplate(config{ipatient}); 
    
    
    % now get some single event data (without overwriting)
    config{ipatient}.LFP.name                                   = {'template1', 'template2', 'template3', 'template4', 'template5', 'template6'};
    LFP{ipatient}                                               = readLFP(config{ipatient}, MuseStruct{ipatient}, false);
    
    ipart = 1;   
    
    % rereference to bipolar
    for markername = ["template1", "template2", "template3", "template4", "template5", "template6"]
        labels_nonum    = regexprep(LFP{ipatient}{ipart}.(markername).label, '[0-9_]', '');
        [~,~,indx]      = unique(labels_nonum);
        for i = 1 : max(indx)
            cfgtemp             = [];
            cfgtemp.reref       = 'yes';
            cfgtemp.refmethod   = 'bipolar';
            cfgtemp.channel     = LFP{ipatient}{ipart}.(markername).label(indx==i);
            group{i}            = ft_preprocessing(cfgtemp,LFP{ipatient}{ipart}.(markername));
        end
        LFP{ipatient}{ipart}.(markername) = ft_appenddata([],group{:});
        clear group
    end
end
    
%     % demean
%     for markername = ["template1", "template2", "'template3", "template4", "template5", "template6"]
%         try
%         cfgtemp             = [];
%         cfgtemp.demean      = 'yes';
%         LFP{ipatient}{ipart}.(markername) = ft_preprocessing(cfgtemp,  LFP{ipatient}{ipart}.(markername));
%         catch
%         end
%     end
    
    
    figure;
    
    maxrange = 0;
    for itemplate = 1 : 6
        temp = vertcat(LFP{ipatient}{ipart}.(sprintf('template%d', itemplate)).trial{:});
        maxrange = max(max(max(abs(temp))), maxrange);
    end
    maxrange = maxrange / 2;
    
    for itemplate = 1 : 6
        subplot(2, 6, itemplate);
        hold;
        indx = randperm(size(LFP{ipatient}{ipart}.(sprintf('template%d', itemplate)).trial, 2), min(100, size(LFP{ipatient}{ipart}.(sprintf('template%d', itemplate)).trial, 2)));
        

        for itrial = indx
            n = 1; ytick = []; label = [];
            for ichan = 1 : size(LFP{ipatient}{ipart}.(sprintf('template%d', itemplate)).label, 1)
                ytick = [ytick, n*maxrange];
                x       = LFP{ipatient}{ipart}.(sprintf('template%d', itemplate)).time{itrial};
                y       = LFP{ipatient}{ipart}.(sprintf('template%d', itemplate)).trial{itrial}(ichan, :);
                lh = plot(x, y + n*maxrange, 'k');
                lh.Color = [lh.Color 0.1];
                label{ichan} = LFP{ipatient}{ipart}.(sprintf('template%d', itemplate)).label{ichan};
                n = n + 1;
            end
        end
        n = 1;
        for ichan = 1 : size(LFP{ipatient}{ipart}.(sprintf('template%d', itemplate)).label, 1)
            x       = LFP_cluster{ipatient}{itemplate}.time;
            y       = LFP_cluster{ipatient}{itemplate}.avg(ichan, :);
            plot(x, y + n*maxrange, 'r', 'linewidth', 2);
            n = n + 1;
        end
        xlim([LFP_cluster{ipatient}{itemplate}.time(1), LFP_cluster{ipatient}{itemplate}.time(end)]);
        ylim([n/2, maxrange * n + n/2]);
    end
    
    
   
end