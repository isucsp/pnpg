
system(sprintf('mv %s  %s.bak',filename,filename));
try
    save(filename);
catch
    keyboard
end

