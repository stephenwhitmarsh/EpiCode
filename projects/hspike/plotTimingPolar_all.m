function plotTimingPolar_all(config, MuseStruct)

t = table;
for ipatient = 1 : size(MuseStruct, 2)
    config{ipatient}.plot.label = {'template1', 'template2', 'template3', 'template4', 'template5', 'template6'};
    
    for ipart = 1 : size(MuseStruct{ipatient}, 2)
        for idir = 1 : size(MuseStruct{ipatient}{ipart}, 2)
            for marker = string(config{ipatient}.plot.label)
                if isfield(MuseStruct{ipatient}{ipart}{idir}.markers, marker)
                    h = size(MuseStruct{ipatient}{ipart}{idir}.markers.(marker), 2);
                    if h > 0
                        t_temp          = table;
                        t_temp.clock    = MuseStruct{ipatient}{ipart}{idir}.markers.(marker).clock';
                        t_temp.marker   = repmat(marker, height(t_temp), 1);
                        t_temp.part     = repmat(ipart, height(t_temp), 1);
                        t_temp.patient  = repmat(ipatient, height(t_temp), 1);
                        t               = [t; t_temp];
                    end
                end
            end
        end
    end
    t_all{ipatient} = t;
end

t.minute = hour(t.clock)*60 + minute(t.clock);
t.theta  = t.minute / (24*60) * 2 * pi;

% save data to table for R
fname   = fullfile(config{ipatient}.datasavedir, 'hypdata_table');
writetable(t, fname);



fig = figure;
subplot(2,6,1);
nbins = 24;
for ipatient = 1 : size(MuseStruct, 2)
    polarhistogram((hour(t_all{ipatient}.clock)*60 + minute(t_all{ipatient}.clock)) / (24*60) * pi*2, nbins*2, 'normalization', 'probability', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    hold on
    title('All days, all templates'); 
    Ax                      = gca;
    Ax.ThetaDir             = 'clockwise';
    Ax.ThetaZeroLocation    = 'top';
    Ax.RTick                = 0;
    Ax.ThetaTick            = linspace(0, 360, nbins+1);
    %     Ax.Alphamap =Ax.Alphamap / 10;
    labels                  = linspace(0, 24, nbins+1);
    thetaticklabels(labels)
end

subplot(2,6,2);
t_sel = t(t.part <= 3, :);
polarhistogram((hour(t_sel.clock)*60 + minute(t_sel.clock)) / (24*60) * pi*2, nbins*2, 'normalization', 'probability');
title('All days, all templates');
Ax                      = gca;
Ax.ThetaDir             = 'clockwise';
Ax.ThetaZeroLocation    = 'top';
Ax.RTick                = 0;
Ax.ThetaTick            = linspace(0, 360, nbins+1);
labels                  = linspace(0, 24, nbins+1);
thetaticklabels(labels)

iplot = 1;
for marker = string(config.plot.label)
    subplot(2,6,6+iplot);

    t_sel = t(t.marker == marker, :);
    polarhistogram((hour(t_sel.clock)*60 + minute(t_sel.clock)) / (24*60) * pi*2, nbins*2, 'normalization', 'probability');
    Ax                      = gca;
    Ax.ThetaDir             = 'clockwise';
    Ax.ThetaZeroLocation    = 'top';
    Ax.RTick                = 0;
    Ax.ThetaTick            = linspace(0, 360, nbins+1);
    labels                  = linspace(0, 24, nbins+1);
    thetaticklabels(labels)

    title(marker);
    iplot = iplot + 1;
end

fname_fig = fullfile(config.imagesavedir, 'polarplots', [config.prefix, '_IED_polar.png']);
isdir_or_mkdir(fileparts(fname_fig));
exportgraphics(fig, fname_fig);
