function [shifted, nshift] = alignXcorr(input,maxiter)

if size(input,1) == 1
    fprintf('Cannot align only one trial!\n');
    return
end

shifted = input;
nshift  = zeros(1, size(input,1));
toshift = zeros(1, size(input,1));

for iter = 1 : maxiter
    avg = nanmean(shifted);
    for itrial = 1 : size(input,1)
        [X2,lags]           = nanxcorr(avg, shifted(itrial,:));
        [~, lagindx]        = findpeaks(X2);
        [~, shiftindx]      = min(abs(lags(lagindx))) ;
        toshift(itrial)     = lags(lagindx(shiftindx));
        nshift(itrial)      = nshift(itrial) + toshift(itrial);
        shifted(itrial,:)   = shift(shifted(itrial,:), toshift(itrial));
    end
    fprintf('Iteration %d of max %d: max absolute shift: %d samples\n', iter, maxiter, max(abs(toshift)));
    if max(abs(toshift)) <= 1
        fprintf('No need to continue\n');
        return
    end
end

