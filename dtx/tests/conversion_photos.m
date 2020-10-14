imagepath = 'Z:\analyses\lgi1\DCIM\100CANON';
imglist = dir2(imagepath);

i_count = 0;
for i_img = 1:size(imglist,1)
    %select only img with the good extension
    if strcmp(imglist(i_img).name(end-3:end), '.CR2')
        i_count = i_count+1;
        fprintf('Converting image %d \n', i_count);
        img_name = fullfile(imagepath, imglist(i_img).name);
        img = imread(img_name);
        imwrite(img, fullfile([img_name(1:end-4), '.png']));
    end
end
