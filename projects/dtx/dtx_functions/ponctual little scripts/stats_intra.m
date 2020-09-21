intra_values = readtable('Z:\analyses\lgi1\Git-Paul\EpiCode\projects\dtx\ponctual little scripts\Intra_values.csv');
CV2_dtx = intra_values.CV2_dtx(~isnan(intra_values.CV2_dtx));
CV2_ctrl = intra_values.CV2_ctrl(~isnan(intra_values.CV2_ctrl));
Freq_DTX = intra_values.Freq_DTX(~isnan(intra_values.Freq_DTX));
Freq_ctrl = intra_values.Freq_ctrl(~isnan(intra_values.Freq_ctrl));

pval_cv2 = ranksum(CV2_dtx, CV2_ctrl);
pval_freq = ranksum(Freq_DTX, Freq_ctrl);