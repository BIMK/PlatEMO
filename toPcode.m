function toPcode()
clc; format compact

    folder = 'GUI\Modules';
    
    files = what(fullfile(cd,folder));
    files = files.m;
    for i = 1 : length(files)
        pcode(fullfile(cd,folder,files{i}));
    end
end