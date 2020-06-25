function compute_spike_xcorr(slurm_task_id)

if ispc
    addpath \\lexport\iss01.charpier\analyses\lgi1\Git-Paul\fieldtrip;
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\projects\dtx'));
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\external'));
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\shared'));
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\NeuralynxMatlabImportExport_v6.0.0'));
    addpath(genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\SPIKY_2_Dez_2019\SPIKY'))
    cd \\lexport\iss01.charpier\analyses\lgi1\Git-Paul\SPIKY_2_Dez_2019\SPIKY\SPIKY_MEX
elseif isunix
    addpath /network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/fieldtrip/
    addpath(genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/projects/dtx'));
    addpath(genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/external'));
    addpath /network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/shared
    addpath (genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/NeuralynxMatlabImportExport_v6.0.0'));
    addpath(genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/SPIKY_2_Dez_2019/SPIKY'))
    cd /network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/SPIKY_2_Dez_2019/SPIKY/SPIKY_MEX
    
end

ft_defaults

%config
config = dtx_setparams_probe_spikes([]);
%table with manual annotation about each unit
if ispc
    unit_table = readtable('\\lexport\iss01.charpier\analyses\lgi1\DTX-PROBE\classification_units.xlsx');
elseif isunix
    unit_table = readtable('/network/lustre/iss01/charpier/analyses/lgi1/DTX-PROBE/classification_units.xlsx');
end

unit_table = table2struct(unit_table);

ipart=1;


%FIXME : dans script projet : script pour filtrer spiketrials, script xcorr
% en fonction

for irat = slurm_task_id
    
    [MuseStruct]                    = readMuseMarkers(config{irat}, false);
    
    if strcmp(config{irat}.type, 'dtx')
        MuseStruct = dtx_remove_wrong_seizure(config{irat}, MuseStruct,true);
        SpikeTrials = readSpikeTrials_MuseMarkers(config{irat}, [],[], false);
    elseif strcmp(config{irat}.type, 'ctrl')
        SpikeTrials = readSpikeTrials_continuous(config{irat}, [],[], false);
    end
    
    [SpikeTrials, ~]                = removetrials_MuseMarkers([], SpikeTrials, MuseStruct); %remove BAD trials
    
    for ilabel = 1:size(config{irat}.name,2)
        
        for i_group = ["mua", "sua"]
            
            spikedata = SpikeTrials{ipart}{ilabel};
            
            %% select units according to their group
            %do it after filtering and keeping only sua
            rat_idx = strcmp({unit_table.ratID}, config{irat}.prefix(1:end-1))';
            rat_table = unit_table(rat_idx,:);
            for i_unit = 1:size(spikedata.label, 2)
                %look all the rat_table unit IDs, select which one is i_unit
                unit_idx = strcmp(split(sprintf('cluster_%d,', rat_table.clusterID), ','), spikedata.label{i_unit});
                % find group info of the unit
                if sum(unit_idx) == 1
                    spikedata.cluster_group{i_unit}              = rat_table(unit_idx).group;
                end
            end
            
            %replace empty cells
            for i_unit = 1:size(spikedata.label, 2)
                if isempty(spikedata.cluster_group{i_unit})
                    spikedata.cluster_group{i_unit} = 'noise';
                end
            end
            
            %filter spikedata
            if strcmp(i_group, "mua")
                group_sel = ~contains(spikedata.cluster_group, 'noise');
            elseif strcmp(i_group, "sua")
                group_sel = ~contains(spikedata.cluster_group, 'noise') & ~contains(spikedata.cluster_group, 'mua');
            end
            group_sel = find(group_sel);
            
            cfgtemp = [];
            cfgtemp.spikechannel = group_sel;
            spikedata = ft_spike_select(cfgtemp, spikedata);
            
            %% compute xcorr
            cfgtemp             = [];
            cfgtemp.binsize     = 0.0005;
            cfgtemp.maxlag      = 0.015;
            cfgtemp.debias      = 'yes';
            cfgtemp.method      = 'xcorr';
            cfgtemp.outputunit  = 'raw';
            xcorr{ilabel}.(i_group)               = ft_spike_xcorr(cfgtemp,spikedata);
            xcorr{ilabel}.(i_group).analysis_name = config{irat}.name{ilabel};
            xcorr{ilabel}.(i_group).group         = i_group;
            %xcorr{ilabel}.(i_group).xcorr : unit * unit * y-values
            
            %% plot xcorr
            fig = figure;
            sgtitle('Legend : cluster_nr(electrode_nr). Top = y axis. Left = x axis.','Interpreter','none','Color','r');
            set(fig, 'units','normalized','position', [0 0 1 0.5]);
            i = 0;
            
            ft_progress('init', 'text',     'Please wait...');
            nr_xcorrs = size(xcorr{ilabel}.(i_group).xcorr,1) * size(xcorr{ilabel}.(i_group).xcorr,2);
            for ix = 1 : size(xcorr{ilabel}.(i_group).xcorr,1) %the unit index for all the line
                for iy = 1 : size(xcorr{ilabel}.(i_group).xcorr,2) %the unit index for all the column
                    
                    i = i + 1; %count subplot
                    
                    ft_progress(i/nr_xcorrs, 'Plotting xcorr %d from %d (rat %d label %d)', i,nr_xcorrs, irat, ilabel);
                    
                    xcorr{ilabel}.(i_group).maxchan(ix) = spikedata.template_maxchan(ix);
                    
                    if ix > iy
                        c = [0 0 0];
                    end
                    
                    if ix < iy
                        c = [0 0 0];
                    end
                    
                    if ix == iy
                        c = [0 0 0.8];
                    end
                    
                    x =xcorr{ilabel}.(i_group).time; %abscisse xcorr
                    y = squeeze(xcorr{ilabel}.(i_group).xcorr(ix,iy,:)); %value of the xcorr : cluster ix in x-axis, iy in y-axis
                    
                    
                    
                    h = subplot(size(xcorr{ilabel}.(i_group).xcorr,1),size(xcorr{ilabel}.(i_group).xcorr,2),i);
                    hold;
                    
                    if ix == 1
                        unit_label = split(xcorr{ilabel}.(i_group).label{iy}, '_');
                        title(sprintf('%s(%d)', unit_label{2},spikedata.template_maxchan(iy)),'FontSize',5);
                    end
                    if iy ==1
                        unit_label = split(xcorr{ilabel}.(i_group).label{ix}, '_');
                        ylabel(sprintf('%s(%d)', unit_label{2},spikedata.template_maxchan(ix)),'FontSize',5);
                        set(gca,'FontWeight','bold');
                    end
                    xticks([])
                    yticks([]);
                    
                    if ~any(isnan(y))
                        
                        %                 Lx = 1:length(x)/2;
                        %                 Rx = length(x)/2 : length(x);
                        %
                        %                 xintL = linspace(x(Lx(1)),x(Lx(end)),100)';
                        %                 yintL = spline(x(Lx),y(Lx),xintL);
                        %                 yintL = smooth1q(yintL,10);
                        %
                        %                 xintR = linspace(x(Rx(1)),x(Rx(end)),100)';
                        %                 yintR = spline(x(Rx),y(Rx),xintR);
                        %                 yintR = smooth1q(yintR,10);
                        
                        
                        bar(x,y, 'FaceColor', c, 'EdgeColor', c);
                        %                 plot(xintL,yintL,'r','linewidth',1);
                        %                 plot(xintR,yintR,'r','linewidth',1);
                        axis tight
                        ax = axis;
                        ylim([0,ax(4)]);
                        pbaspect([1 1 1])
                        
                        
                        %                 axis off
                    end
                end
            end
            ft_progress('close')
            
            %% print to file
            fprintf('Print to file and save data\n');
            
            if ~isfolder(fullfile(config{irat}.imagesavedir,'xcorrs'))
                mkdir(fullfile(config{irat}.imagesavedir,'xcorrs'));
            end
            
            fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
            set(fig,'PaperOrientation','landscape');
            set(fig,'PaperUnits','normalized');
            set(fig,'PaperPosition', [0 0 1 1]);
            print(fig, '-dpdf', fullfile(config{irat}.imagesavedir,[config{irat}.prefix,'p',num2str(ipart),'-xcorr_',str2mat(i_group),'_',config{irat}.name{ilabel},'.pdf']),'-r600');
            print(fig, '-dpng', fullfile(config{irat}.imagesavedir,[config{irat}.prefix,'p',num2str(ipart),'-xcorr_',str2mat(i_group),'_',config{irat}.name{ilabel},'.png']),'-r600');
            
            close all
            
        end %i_group
    end %ilabel
    
    fname = fullfile(config{irat}.datasavedir, [config{irat}.prefix(1:end-1), '_xcorr.mat']);
    save(fname, 'xcorr', '-v7.3');
    clear xcorr
end %irat

end