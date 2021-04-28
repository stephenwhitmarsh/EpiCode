function out = loadload(fname)

count = 0;
err_count = 0;
while count == err_count
    try
        out = load(fname);
    catch ME
        err_count = err_count + 1;
    end
    count = count + 1;
end