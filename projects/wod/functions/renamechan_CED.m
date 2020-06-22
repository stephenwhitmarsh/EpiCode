function chantitle  = CEDrename_chan(chantitle,channr,otherchanstitles,doprint)
%DESCRIBE.
% if has already the name : first chan keep the name
% otherchanstitles can be []

% cant make fieldnames with minusses
if any(ismember('-',chantitle))
    chantitle = strrep(chantitle,'-','_');
    if doprint
        fprintf('Channel %s is renamed %s (cant make fieldnames with minusses)\n', chantitle, strrep(chantitle,'-','_'));
    end
end

% cant make fieldnames with white spaces
if any(ismember(' ',chantitle))
    chantitle = strrep(chantitle,' ','_');
    if doprint
        fprintf('Channel %s is renamed %s (cant make fieldnames with white spaces)\n', chantitle, strrep(chantitle,' ','_'));
    end
end

%create name if is empty
if isempty(chantitle)
    chantitle = sprintf('chan%d', channr);
    if doprint
        fprintf('Channel %d is renamed %s (has no name)\n', channr, chantitle);
    end
end

%add 'x' to the field name if it begins with a number
if ismember(chantitle(1), ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9'])
    chantitle = insertBefore(chantitle,chantitle,'X');
    if doprint
        fprintf('Channel %s is renamed %s (begins with a number)\n', chantitle(2:end), chantitle);
    end
end

%rename channel if it has the same name that a previous one
if any(strcmp(chantitle, otherchanstitles))
    oldname = chantitle;
    chantitle = sprintf('%s_chan%d', chantitle, channr);
    if doprint
        fprintf('Channel %d %s as the same name than a previous channel : renamed %s\n', channr, oldname, chantitle);
    end
end

end

