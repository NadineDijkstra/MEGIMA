function [session, sessionName] = current_session(filename, session, dataPath)

%[session, sessionName] = current_session(filename, session, dataPath)
%
% Returns the number and name of the current session. 
% If session name already exists, user will be prompted for action.


session = sprintf('%02d', session);
sessionName = [filename '_' session];

if exist(fullfile(dataPath, sprintf('%s.mat',sessionName)), 'file')
    query = input('\nWARNING: Output filename already exists. \nDo you want to specify a different run number? (y or n) ','s');
	if query(1) == 'n'
        disp(['File will be temporarily saved as "temp' sessionName '.mat" Please rename the file!']);
        sessionName = ['temp_' sessionName];
    elseif query(1) == 'y'
        session = input('What run number is this?', 's');
        sessionName = [filename '_' sprintf('%02d', str2double(session))];
    else
        error('No valid input given');
	end
end

session = str2double(session);
disp(['Running session: ' sessionName]); disp(' ');

