figure;hold;
pn_mua = scatter(1,1,'^k');
pn_sua = scatter(1,1,'^k', 'filled');
in_mua = scatter(1,1,'ok');
in_sua = scatter(1,1,'ok', 'filled');
legend([pn_mua, pn_sua, in_mua, in_sua], 'PN (MUA)', 'PN (SUA)', 'IN (MUA)', 'IN (SUA)');