function [seqevents, epoch_times] = read_events4(infile); 
% READ_EVENTS4  reads Events.txt file from cheetah (via dat2txt.exe conversion)
%              interprets port values in terms of sequence experiment
%              events.  Written for BX-35 controlled programs, starting
%              7/04.
%              
% 			[seqevents, epoch_times] = read_events(infile);
%
% infile can be omitted and program will look for events.txt in current directory
% seqevents is an array with the following columns:
%
%        seqevents = [timestamp light feeder stim zone];
%
% timestamps are in microseconds (standard cheetah format)
% light, feeder and stim are binary, 1 = on, 0 = off
% zone numbers range from 0 to 9
% note that tone events show up as [ts 0 0 0 zone]
% 
% epoch_times is an array of the start and stop times of each behavioral epoch
% they are in microseconds (same as cheetah timestamps)
% these epochs are denoted by "a" and "b" events inserted into the events.dat file
% during recording session.  
% epoch_times = 
% [start1 stop1]
% [start2 stop2]
% ...
% 7/26/2004 DE altered to accomodate new experiment with BX-35, includes
% light onset info.
% 5/8/2003  DE fixed so that bits 9 and 10 are also read and translated into light events 8 & 9
% 4/17/2003 DE accelerated by having it read the whole file into a string first
%           this is a revision to read_events2x
% 3/3/2003 DE added ability to read old AND new neuralynx event files


if nargin<1
	infile = 'events.txt';
end;

fid = fopen(infile);
if fid==-1
	error(['Could not open file: ' infile]);
end;
	
end_of_line_char = 10;

i = 1;
epoch_row = 0;
last_flag = '';

is_new = [];  % used to deterime if file is of new cheetah format or old


[newdata readcnt] = fread(fid, inf, 'uchar');  % read the entire file into a string
filedata = char(newdata);
fclose(fid);
line_ends = find(filedata==end_of_line_char);
start_strpos = 1;

% pre-allocate storage for faster execution
num_lines = length(line_ends);
timestamp = zeros(num_lines,1);
light = zeros(num_lines, 1);
feeder = zeros(num_lines,1);
stim = zeros(num_lines,1);
zone = zeros(num_lines,1);


for j = 1:length(line_ends)
   
   curline = filedata(start_strpos:line_ends(j)-2)';  % note we skip back two places to avoid the linespace and carriage return that are both used in dos text files for a return
   %disp(curline)
   start_strpos = line_ends(j)+1;
   
   % deal with header and blank lines
	if length(curline)<5 | all(curline(1:4)=='none')
		continue; 
	end; 
	[ts] = sscanf(curline, '%*i,%*i,%*i,%lf');
	[eventstr] = sscanf(curline, ...
		'%*i,%*i,%*i,%*lf,%*i,%*i,%*i,%*i,%*i,%*i,%*i,%*i,%*i,%*i,%*i,%*i,%*i,%41c');
    [m n] = size(eventstr);
    if m>n  % needed for version matlab 7.5
        eventstr = char(eventstr');
    end;
	ttltxti = findstr(eventstr, 'TTL Value:');
   
	if ~isempty(ttltxti) & length(eventstr)>=38
       is_new = findstr(eventstr(ttltxti+10:length(eventstr)), 'x');
       if ~isempty(is_new)
           ttlvalstr = eventstr(ttltxti+13:length(eventstr));  % note: values look like:
																			 % 0xE080, we strip off the 0x
       else
           ttlvalstr = eventstr(ttltxti+15:length(eventstr));  % note: values look like:
																			 % FFFFE080, we strip off
																			 % the four leading Fs
       end
       
       ttlvalbin = dec2bin(hex2dec(ttlvalstr),16);  % (e.g., 1100000000100000)
		%disp(ttlvalbin);
	   %valbinflip = fliplr(ttlvalbin);              % reverse the binary number
													% so that string indices correspond
										            % to bit locations
      % the first four bits encode the zone number in binary format
      % the 5th bit is light
      % the 6th bit is feeder
      % the 7th bit is stim
      zone(i) = bin2dec(ttlvalbin(13:16));
      %disp([ttlvalbin 'zone: ' num2str(zone(i))])
      if ttlvalbin(17-5)=='1'
         light(i) = 1;
      else
         light(i) = 0;
      end;      
      if ttlvalbin(17-6)=='1'
         feeder(i) = 1;
      else
         feeder(i) = 0;
      end;
      if ttlvalbin(17-7)=='1'
         stim(i) = 1;
      else
         stim(i) = 0;
      end;
      
      timestamp(i) = ts;
	   i = i+1;
   else  % if ~isempty(ttltxti) & length(eventstr)>=38
      % NOW PROCESS ANY NON-TTL INPUT
       % Process the User entered event codes which denote epochs

       cur_event_flag = sscanf(eventstr, '%s');  % remove white space
       switch cur_event_flag 
			case 'a'
              if ~isempty(last_flag) & (last_flag == 'a')
                 disp(['Warning: Missing "b" event flag before "a" at timestamp ' num2str(ts)])
              end;
              epoch_row = epoch_row + 1;
              epoch_times(epoch_row, 1) = ts;
              last_flag = 'a';
			case 'b'
              if last_flag ~= 'a' 
                 disp(['Warning: "b" event flag without "a" at timestamp ' num2str(ts)])
                 epoch_row = epoch_row + 1;
              end;
              epoch_times(epoch_row, 2) = ts;
              last_flag = 'b';
			otherwise
              disp(['Warning: Unknown User Input Event, "' eventstr '" , at timestamp ' num2str(ts)]) 
       end; 
       
    end;  % if ~isempty(ttltxti) & length(eventstr)>=38

end;

%if ~isempty(is_new)
%    disp('File is a new cheetah file')
%else
%    disp('File is an old cheetah file');
%end

% CUT OFF NON-ALLOCATED VALUES FROM THE END OF LISTS
% note that we pre-allocated space for these variables, 
% set to the number of lines in the file.  
% Since not every line is a valid event, there are extra zeros 
% at the end of every variable.  Here, we remove those zeros.
timestamp = timestamp(1:i-1,:);
feeder = feeder(1:i-1, :);
stim = stim(1:i-1, :);
zone = zone(1:i-1, :);
light = light(1:i-1, :);

seqevents = [timestamp light feeder stim zone];

% remove any 'all clear' events (sent to output port to turn off lights)
keepi = ~(feeder==0 & stim==0 & zone==-1);
seqevents = seqevents(keepi,:);
