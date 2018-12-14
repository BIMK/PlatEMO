function search()
clc;format compact

    word  = 'Sciences';
    exist = true;

    Folders = strsplit(genpath(fileparts(mfilename('fullpath'))),';');
    for i = 1 : length(Folders)
        folder = what(Folders{i});
        if ~isempty(folder)
            files = folder.m;
            for j = 1 : length(files)
                file = fullfile(Folders{i},files{j});
                f    = fopen(file);
                line = 0;
                Exi  = false;
                while ~feof(f) && ~Exi
                    line = line + 1;
                    if ~isempty(regexp(fgetl(f),['\<',word,'\>'],'once'))
                        Exi = true;
                    end
                end
                if exist == false && Exi == false
                    fprintf('<%s> 不存在字符串 %s\n',file,word);
                elseif exist == true && Exi == true
                    fprintf('<%s> 的 <第%d行> 存在字符串 %s\n',file,line,word);
                end
                fclose(f);
            end
        end
    end
end