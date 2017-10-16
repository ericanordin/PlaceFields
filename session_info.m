function session_info_struct = session_info(rat, session_num, depth_flag);
% SESSION_INFO  returns critical info about each data session
%        (reads from session_data.xls file)
%   info = session_info(rat_num, session_num, depth_flag=0)
%          if depth flag is non-zero, depths will be extracted, otherwise
%          they will not.  (depths take a while)
%   info is a structure with the following fields:
%           sess_num
%           data_path
%           tfile_path
%           heart_file
%           events_file
%           events_reader  this is the name of the function that will parse
%                          the events file
%           pos_file1     this is the position file (e.g., vtm1.pvd)
%           pos_file2
%           pos_sleep_file1   position during sleep (e.g., vts1.pvd)
%           pos_sleep_file2
%           pos_sleep_file3
%           pos_ts_flag   if non-zero, indicates pvd tracker file used ts units (100 usec) not usec 
%           no_maze       number of maze sessions 
%           ini_file1
%           ini_file2
%           type          'V', 'S', 'V2S', 'S2V', 'V2V', 'prob' --these are
%                          describe the type of experiment run that day
%           days           days experience with the particular sequence
%           tet_depth     array of tetrode depths in microns
%
% Requires the excel file session_data.xls to be in matlab path.  See the
% file for headers which tell what each column contains.
% 
% Here are more details on the type and days fields:
%
% V   = standard complex sequence with two repeated segments
% prob= simple sequence with probablistic reward and extra zones
% S   = standard simple sequence with 6 elements
% S2V = transition from simple to complex V
% V2S = transition from complex V to simple
% V2V = one day that i transfered complex sequence to another complex sequence
% TW = T-maze Woods Task
% TC = T-maze Continuous (rewards only at ends of T arms)
% The number indicates days of experience with that particular sequence
%
% note that I count the first transition session as day 1 for the second 
% sequence.  Thus, V2S counts as day one for the simple sequence (and I
% ignore the V part because its only half of the first maze session).
% For prob sessions, I count prior experience with the same underlying simple
% sequence, even though the prob extra is technically a different sequence
%
% Note also that this function uses persistent variables to speed up
% execution time.  If you edit the exel spreadsheet 'session_data.xls' you
% will need to clear the function in order to cause reload of file.  Just
% type:
%    clear session_info

% this is info on the type of experiment run (see above for more details)
% its format is [session_num   session_type   days_training_on_same_sequence]

persistent sess_data_num sess_data_txt

if nargin<=2
   depth_flag = 0;
end;

if isempty(sess_data_txt)
    disp('reading session_data.xls file');
    [sess_data_num sess_data_txt] = xlsread('session_data.xls'); 
    sess_data_txt = sess_data_txt(2:end, :);  % delete first line, which is the header.  
    %sess_data_num = sess_data_num(2:end, :); 
end;

% SET SESSION AND RAT NUM
session_info_struct.sess_num = session_num;
session_info_struct.rat = rat;

% DETERMINE PATH
folder_name = [num2str(rat) '_' sprintf('%02d', session_num) 'p'];
%folder_name = [num2str(rat) '_' sprintf('%02d', session_num)];
folder_found = 0;
if exist([folder_name], 'dir')
   data_path_base = pwd;
   folder_found = 1;
end;
if exist(['d:\data\' folder_name], 'dir')
   data_path_base = 'd:\data';
   folder_found = 1;
end;
if exist(['f:\data\' folder_name], 'dir')
   data_path_base = 'f:\data';
   folder_found = 1;
end;
if ~folder_found & exist(['f:\raw_data\' folder_name], 'dir')
   data_path_base = 'f:\raw_data';
   folder_found = 1;
end;  
if ~folder_found & exist(['d:\vidanalysis\' folder_name], 'dir')
   data_path_base = 'd:\vidanalysis';
   folder_found = 1;
end;  
% if ~folder_found & exist(['\\euston-hpxw8400\Data\' folder_name], 'dir')
%    data_path_base = '\\euston-hpxw8400\Data';
%    folder_found = 1;
% end; 

% if ~folder_found & exist(['c:\data\vidanalysis\' folder_name], 'dir');
%    data_path_base = 'c:\data\vidanalysis';
%    folder_found = 1;   
% end;
% if exist(['e:\data\' folder_name], 'dir')
%    data_path_base = 'e:\data';
%    folder_found = 1;
% end;
% if ~folder_found & exist(['e:\data\vidanalysis\' folder_name], 'dir');
%    data_path_base = 'e:\data\vidanalysis';
%    folder_found = 1;   
% end;
if ~folder_found & exist(['t:\prefrontal\data\' folder_name], 'dir')
   data_path_base = 't:\prefrontal\data';
   folder_found = 1;
end;   
if ~folder_found
   disp(['Warning: Can not find folder: ' folder_name])
   return;
end;
   
%host = hostid;
%if str2num(host{1}) == 153609  % this sets data path local on BBClust6
%   data_path_base = 'c:\data';
%end

data_path = [data_path_base '\' folder_name];

% DESCRIBE TYPE OF EXPERIMENT, DAYS RUN AND NUMBER OF MAZE SESSIONS
sess_i = find(sess_data_num(:,1) == rat & sess_data_num(:, 2) == session_num);
if isempty(sess_i)
    warning(['Session type not found in database for rat ' num2str(rat) ' session ' num2str(session_num)]);
    return;
end;
session_info_struct.type = sess_data_txt{sess_i, 3};
session_info_struct.days = sess_data_num(sess_i, 4);
session_info_struct.no_maze = sess_data_num(sess_i, 5);
session_info_struct.pos_ts_flag = sess_data_num(sess_i, 8);

% CREATE STANDARD FILE PATHS
session_info_struct.data_path = data_path;
session_info_struct.tfile_path = [data_path '\tfiles'];

if exist([data_path '\heartbeat.t'], 'file')
   session_info_struct.heart_file = [data_path '\heartbeat.t'];
else
   session_info_struct.heart_file = '';
end;
session_info_struct.events_file = [data_path '\events.txt'];

% FIND RAW EVENT FILE (USED LATER BECAUSE IT CONTAINS THE DATE NEEDED FOR
% DEPTH INFO)
if depth_flag
    depth_ok = true;
   if exist([data_path '\events.dat'], 'file')
       raw_ev_file = [data_path '\events.dat'];
	elseif exist([data_path '\events.nev'], 'file')
       raw_ev_file = [data_path '\events.nev'];
    elseif exist([data_path '\EventRecords.dat'], 'file')
       raw_ev_file = [data_path '\EventRecords.dat'];
       depth_ok = false;  % cannot read date from this file and not depth file exists anyway
	else
       disp(['Warning: Raw events file not found for rat ' num2str(rat) ' session ' num2str(session_num)]);
       disp(['Depth information not available']);
       depth_ok = false;
	end;
end;

% HANDLE VIDEO FILE INFO, INCLUDING A SEARCH FOR MAZE 2 VIDEO
if rat==6862 & session_info_struct.pos_ts_flag==1
   session_info_struct.pos_file1 = [data_path '\vt.pvd'];
   session_info_struct_pos_file2 = '';
else  % for all other rats
	session_info_struct.pos_file1 = [data_path '\vtm1.pvd'];
	if ~exist([data_path '\vtm1.pvd'], 'file')
        disp(['Video file does not exist: ' session_info_struct.pos_file1])
	end;
	if exist([data_path '\vtm2.pvd'], 'file')
        session_info_struct.pos_file2 = [data_path '\vtm2.pvd'];
	else
        session_info_struct.pos_file2 = '';
	end;
end;

if exist([data_path '\vts1.pvd'], 'file')
     session_info_struct.pos_sleep_file1 = [data_path '\vts1.pvd'];
else
     session_info_struct.pos_sleep_file1 = '';
end;
if exist([data_path '\vts2.pvd'], 'file')
     session_info_struct.pos_sleep_file2 = [data_path '\vts2.pvd'];
else
     session_info_struct.pos_sleep_file2 = '';
end;
if exist([data_path '\vts3.pvd'], 'file')
     session_info_struct.pos_sleep_file3 = [data_path '\vts3.pvd'];
else
     session_info_struct.pos_sleep_file3 = '';
end;

% SPECIFY INI FILES   
session_info_struct.ini_file1 = [data_path '\' sess_data_txt{sess_i,6} ];
if ~isempty(sess_data_txt{sess_i,7})
    session_info_struct.ini_file2 = [data_path '\' sess_data_txt{sess_i,7}];
else
    session_info_struct.ini_file2 = '';
end;

% SPECIFY EVENT FILE READER FUNCTION
session_info_struct.event_reader = sess_data_txt{sess_i, 9};

if depth_flag

    if depth_ok
        % GET DATE FROM EVENT FILE
        [header] = textread(raw_ev_file, '%[^\n]', 3);  % read first three lines from file
        date_str = '';
        for i = 1:3
           start_i = findstr(header{i}, 'Time Opened: (m/d/y):');
           if ~isempty(start_i)
              date_str = char(strread(header{i}(start_i+21:end), '%s', 1));
              break;
           end
        end;
	
        % GET DEPTHS
        if isempty(date_str)     
           warning('Date cant be extracted from raw Events file');
           session_info_struct.tet_depths = [];
        else
           session_info_struct.tet_depths = read_depth([data_path_base '\' num2str(rat) '_depths.txt'], date_str);
        end;
        
    else
           session_info_struct.tet_depths = [];
    end;
end