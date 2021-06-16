function saveMarker_TFR(data, markername, fname_out)

TFR{size(data, 2)} = [];
for ipart = 1 : size(data, 2)
    try
        TFR{ipart}.(markername) = data{ipart}.(markername);
    catch
    end
end

save(fname_out, 'TFR', '-v7.3')
