function datout = shift(dat,steps)

datout = circshift(dat,steps);

if steps > 0
    datout(1:steps) = nan;
elseif steps < 0
    datout(end+steps:end) = nan;
end

