function [paths,names] = filesearch_substring(startDir,searchstring,norecurse)
    if (~exist('norecurse','var') || isempty(norecurse))
        norecurse=0;
    end
    [paths names] = dir_search(startDir,searchstring,norecurse);
    
    function [paths names] = dir_search(currDir,searchstring,norecurse)
        if nargin < 2
            fprintf('Usage: [paths names] = dirsearch(currDir, searchstring)\n');
            paths = '';
            names = '';
            return;
        elseif (exist('norecurse','var')~=1)
            norecurse=0;
        end
        paths = {};
        names = {};
        list_currDir = dir(currDir);

        for u = 1:length(list_currDir)
            if (list_currDir(u).isdir==1 && strcmp(list_currDir(u).name,'.')~=1 && strcmp(list_currDir(u).name,'..')~=1 && norecurse==0)
                [temppaths tempnames] = dir_search(sprintf('%s%s%s',currDir,filesep,list_currDir(u).name),searchstring);
                paths = {paths{:} temppaths{:}};
                names = {names{:} tempnames{:}};
            elseif (length(list_currDir(u).name) > 4)
                    if isempty(strfind(list_currDir(u).name, searchstring))==0
                        paths = {paths{:} currDir};
                        names = {names{:} list_currDir(u).name};
                    end
            end
        end

	end
end

