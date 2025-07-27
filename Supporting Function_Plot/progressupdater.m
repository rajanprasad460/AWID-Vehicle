function str = progressupdater(i,row,str_maj)
% function str = progressupdater(i,row,str_maj)
% progressupdater reads the current loop count, total loop, and user
% defined string to be displayed as an update to the user

s_progress=blanks(50);
progress=round(i/row*50);
progress2 = (i/row*50);
s_progress(1:progress)=repelem('*',progress);
a_progress=sprintf('%.2f%%',progress2*2);
s_progress(21:21+length(a_progress)-1)=a_progress;
%------------------------------------------------------
str{1}=str_maj;
str{2}=sprintf('[%s] Complete \n',s_progress);
% This part clears the displayed so that it looks similar to carriage
% return
if i>1
    %[warnMsg, warnId] = lastwarn;
    [~, warnId] = lastwarn;
    if ~isempty(warnId)
        % Don't delete the characters printed
        lastwarn('') % Clear last warning message
    else % Delete the cahracters printed
        % Deleting the data Written using this function
        fprintf(repmat('\b',1,length(str{1})+length(str{2})+1));
    end
end

%------ Printing the result -------------------------------------
fprintf('%s\n%s',str{1},str{2});


end