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

config = dtx_setparams_probe_spikes([]);

ipart=1;

for irat = 1:7
    fname = fullfile(config{irat}.datasavedir, [config{irat}.prefix(1:end-1), '_xcorr.mat']);
    stats{irat} = load(fname);
end

% good xcorr for each rat
good_xcorr{1} = {[16, 61],[16, 125],[16, 127],[16, 129],[47, 60],[47, 61],[47, 79],[47, 87]}; %noise 64, 124, 51, 82, 84, 86, 89, 96
good_xcorr{2} = {};%noise 8
good_xcorr{3} = {[3, 82],[3, 85],[9, 83],[44, 83],[44, 86]}; %noise 45
good_xcorr{4} = {[4, 24],[4, 19],[4, 12],[4, 7],[7, 23],[7, 24],[8, 26],[9, 43],[16, 24]};%noise 25
good_xcorr{5} = {[0, 20],[1, 9],[1, 17],[1, 21],[1, 25],[1, 30],[1, 76],[9, 25],[11, 78],[14, 77],[15, 74],[15, 77]};%noise 4, 57
good_xcorr{6} = {[15, 16],[15, 17]};%noise 3,4
good_xcorr{7} = {[20, 24],[23, 30],[27, 30],[30, 70]};%noise 1,5,11


for irat = 1:7
  
    if isempty(good_xcorr{irat})
        continue
    end
    
    %re-find cluster labels
    cluster_labels = squeeze(split(stats{irat}.xcorr.label, '_'))';
    cluster_labels = cellfun(@str2num,cluster_labels(2,:));
    
    for i_xcorr = 1:size(good_xcorr{irat},2)
        for i_flip_xcorr = 1:2
            
            %plot x versur y, and y versur x
            if i_flip_xcorr == 2
                good_xcorr{irat}{i_xcorr} = flip(good_xcorr{irat}{i_xcorr});
            end
            
            
            idx_cluster_x = find(good_xcorr{irat}{i_xcorr}(1) == cluster_labels);
            idx_cluster_y = find(good_xcorr{irat}{i_xcorr}(2) == cluster_labels);
            
            if length(idx_cluster_x) > 1 || length(idx_cluster_y) > 1
                error('it should have only one label per unit');
            end
            
            fig = figure;
            
            
            
            x = stats{irat}.xcorr.time; %abscisse xcorr
            y = squeeze(stats{irat}.xcorr.xcorr(idx_cluster_x,idx_cluster_y,:)); %value of the xcorr : cluster ix in x-axis, iy in y-axis
            
            bar(x,y, 'FaceColor', 'k', 'EdgeColor', 'k');
            axis tight
            ax = axis;
            ylim([0,ax(4)]);
            pbaspect([1 1 1])
            
            xlabel(sprintf('%s (chan %d)',stats{irat}.xcorr.label{idx_cluster_x}, stats{irat}.xcorr.maxchan(idx_cluster_x)), 'Interpreter','none');
            ylabel(sprintf('%s (chan %d)',stats{irat}.xcorr.label{idx_cluster_y}, stats{irat}.xcorr.maxchan(idx_cluster_y)), 'Interpreter','none');
            
            if ~isfolder(fullfile(config{irat}.imagesavedir,'xcorrs_good'))
                mkdir(fullfile(config{irat}.imagesavedir,'xcorrs_good'));
            end
            
            fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
            set(fig,'PaperOrientation','landscape');
            set(fig,'PaperUnits','normalized');
            set(fig,'PaperPosition', [0 0 1 1]);
            print(fig, '-dpdf', fullfile(config{irat}.imagesavedir,'xcorrs_good',[config{irat}.prefix,'xcorr_x',num2str(good_xcorr{irat}{i_xcorr}(1)),'_y',num2str(good_xcorr{irat}{i_xcorr}(2)),'.pdf']),'-r600');
            print(fig, '-dpng', fullfile(config{irat}.imagesavedir,'xcorrs_good',[config{irat}.prefix,'xcorr_x', num2str(good_xcorr{irat}{i_xcorr}(1)),'_y',num2str(good_xcorr{irat}{i_xcorr}(2)),'.png']),'-r600');
            
            close all
        end
    end
    
    
end %irat

end