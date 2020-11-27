function saveMarker(data, markername, fname_out)

LFP{size(data, 2)} = [];
for ipart = 1 : size(data, 2)
    try
        LFP{ipart}.(markername) = data{ipart}.(markername);
    catch
    end
end

save(fname_out, 'LFP', '-v7.3')
