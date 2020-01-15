function [shifted,nshift] = align_data(input,maxiter)

shifted = input;
nshift  = zeros(1,size(input,1));

% figure; hold;
for iter = 1 : maxiter
    fprintf('Iteration: %d\n',iter);
    if size(input,1) > 1
    avg = nanmean(shifted);
%     plot(avg);
    
    for itrial = 1 : size(input,1)
%         avg = nanmean(shifted);  % much slower and no visible improvement
        [X2,lags]           = nanxcorr(avg,shifted(itrial,:));
        [~, lagindx]        = findpeaks(X2);
        [~, shiftindx]      = min(abs(lags(lagindx))) ;
        toshift             = lags(lagindx(shiftindx));
        nshift(itrial)      = nshift(itrial) + toshift;
        shifted(itrial,:)   = shift(shifted(itrial,:),toshift);
    end
    else
        fprintf('Oops, cannot align only one trial!\n');
    end
end