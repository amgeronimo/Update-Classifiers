function [files,viability,errorcode, sessions,runs,dates] = FindFiles(handle)

%Could find a way to remove the loops in the initial code, but it shouldnt
%cause too much of a problem time-wise.


dataloc = [handle.dataloc '\data'];
code = handle.participant_id;
allfiles = dir(dataloc);
files = {};
sessions = {};
runs = {};
dates = {};
userlocs = zeros(length(allfiles),1);
waitboxtext = 'Finding Files';
set(handle.listbox1,'String',waitboxtext);
for i = 1:length(allfiles)
    userlocs(i) = strcmp(allfiles(i).name(1:end-3),code);
    if userlocs(i) == 1
        sess = allfiles(i).name(end-2:end);
        %         runs{k} = {};
        userfiles = dir([dataloc '\' allfiles(i).name]);
        for j = 1:length(userfiles)
            runlocs = regexp(userfiles(j).name,'[(.\)]+dat');
            %               = strfind(userfiles(j).name,'.dat');
            if ~isempty(runlocs) %current file a data file
                rstart = regexp(userfiles(j).name, 'R+[0-9]+[(.\)]+');
                if isempty(rstart) %bad data not in Rxx. format
                else
                    sessions = [sessions sess];
                    runs = [runs; userfiles(j).name(rstart+1:runlocs-1)];
                    files = [files; {['\' allfiles(i).name '\' userfiles(j).name]}];
                    Finfo = dir([dataloc '\' allfiles(i).name '\' userfiles(j).name]);
                    dates = [dates; Finfo.date];
                end
            end
        end
    end
end


%Check files for viability
viability = zeros(length(files),1);
errorcode = cell(length(files),1);
modval = length(files)/20;
for i = 1:length(files)
    try
        if handle.ctype == 6
            [viability(i) errorcode{i}] = CheckFileViabilityMI([dataloc files{i}],handle.params);
        else
            [viability(i) errorcode{i}] = CheckFileViabilityP300_v2([dataloc files{i}],handle.params);
        end
    catch
        disp(['error in file ' num2str(i) ', cannot process']);
    end
    try
        if mod(i,ceil(modval))==0
            waitboxtext = [waitboxtext '.'];
            set(handle.listbox1,'String',waitboxtext);
        end
    catch
        disp('error in dotting?');
    end
    pause(.01);
    
end



    
end


