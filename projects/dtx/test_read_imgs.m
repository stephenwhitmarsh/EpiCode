fig = figure;
%subplot (4,2,1)
imagepath = '\\lexport\iss01.charpier\analyses\lgi1\DTX-PROBE\Analyses_Paul\DTX5\spikerate-thresh12-merge_after_tempmatch\spikerate-thresh12-merge_after_tempmatch\DTX5-part1-spikerates_temp_24_pattern_InterIctal.png';
image(imread(imagepath));
test = imread(imagepath);
%hauteur x largeur x 3 couleurs
image(test(3300:4900,4300:6100,:));
axis off

set(fig,'PaperOrientation','landscape');
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition', [0 0 1 1]);
print(fig, '-dpdf', fullfile(cfg.imagesavedir,'test_read_image'),'-append');
print(fig, '-dpng', fullfile(cfg.imagesavedir,'test_read_image'),'-r600');
