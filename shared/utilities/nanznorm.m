function datz = nanznorm(datin)

for rowi = 1 : size(datin,1)
    dat = datin(rowi,:);
    if any(isnan(dat(:)))
        xmu=nanmean(dat);
        xsigma=nanstd(dat);
        datz(rowi,:) = (dat-repmat(xmu,size(dat,1),1))./repmat(xsigma,size(dat,1),1);
    else
        datz(rowi,:) = zscore(dat);
    end
    
end
