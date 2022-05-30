function patch_std(x, ymean, ystd, color)

if nargin < 4
    color = 'k';
end

y = [ymean - ystd; ystd; ystd]';
filled_SD = area(x,y);
filled_SD(1).FaceAlpha = 0; filled_SD(2).FaceAlpha = 0.4; filled_SD(3).FaceAlpha = 0.4;
filled_SD(1).EdgeColor = 'none'; filled_SD(2).EdgeColor = 'none'; filled_SD(3).EdgeColor = 'none';
filled_SD(2).FaceColor = color; filled_SD(3).FaceColor = color;
filled_SD(1).ShowBaseLine = 'off';