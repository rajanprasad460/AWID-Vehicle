function filename = AutoRename(Folder,f_name)
% function filename = AutoRename(Folder,f_name)
% If file exist, create a new file with - [1 2 3 ...] so on
% Folder = cd; % The directory where file has to be stored
% f_name = 'Download_Code_V.bat'; % Filename with extension


f_check = 1;
f_count = 1;
%% Get file parts here
[~,f_name_def,extension] = fileparts(f_name);

%% Check if the file exist and add version numeric
while (f_check)
    if isfile(fullfile(Folder,strcat(f_name_def,extension)))
        f_name_def = strcat(f_name_def,'-',num2str(f_count));
        f_check = 1;
        f_count = f_count+1;
    else
        f_check = 0;
    end
end
%--------------------------------
filename = fullfile(Folder,strcat(f_name_def,extension));