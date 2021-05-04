function patch_std(x, ymean, ystd)

y = [ymean - ystd; ystd; ystd]';
filled_SD = area(x,y);
filled_SD(1).FaceAlpha = 0; filled_SD(2).FaceAlpha = 0.4; filled_SD(3).FaceAlpha = 0.4;
filled_SD(1).EdgeColor = 'none'; filled_SD(2).EdgeColor = 'none'; filled_SD(3).EdgeColor = 'none';
filled_SD(2).FaceColor = 'b'; filled_SD(3).FaceColor = 'b';
filled_SD(1).ShowBaseLine = 'off';