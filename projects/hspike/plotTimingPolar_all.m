function plotTimingPolar_all(config, MuseStruct)

t = table;
for ipatient = 1 : size(MuseStruct, 2)
    config{ipatient}.plot.label = {'template1', 'template2', 'template3', 'template4', 'template5', 'template6'};
    
    for ipart = 1 : size(MuseStruct{ipatient}, 2)
        for idir = 1 : size(MuseStruct{ipatient}{ipart}, 2)
            for marker = string(config{ipatient}.plot.label)
                fprintf('Concatinating Patient %d, part %d, dir %d, marker %s\n', ipatient, ipart, idir, marker);
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
end

t.minute = hour(t.clock)*60 + minute(t.clock);
t.theta  = t.minute / (24*60) * 2 * pi;

% save data to table for R
fname   = fullfile(config{ipatient}.datasavedir, 'hypdata_table');
writetable(t, fname);



fig = figure;
nbins = 24;
polarhistogram((hour(t.clock)*60 + minute(t.clock)) / (24*60) * pi*2, nbins*2, 'normalization', 'probability', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fname_fig = fullfile(config{1}.imagesavedir, 'polarplot_all.png');
isdir_or_mkdir(fileparts(fname_fig));
exportgraphics(fig, fname_fig);
