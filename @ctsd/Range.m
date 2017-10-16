function R = Range(tsa, tsflag)

% ctsd/Range
% 
%  R = Range(tsa, 'sec')
%  R = Range(tsa, 'ts')
%  R = Range(tsa. 'all_ts')
%
%  returns range covered by tsa
%
% INPUTS
%     tsa 
%      tsflag: if 'ts' returns time in timestamps,
%              if 'sec' returns time in sec
%              if 'sec0' returns time in sec counting from 0
%              if 'ms' returns time in ms
%
% ADR 
% version L5.0
% v4.1 28 oct 1998 flag no longer optional
% v5.0 set up to expect inputs with ts in SEC (as per new ADRLAB policy)
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0.
% Version control M3.0.

R = StartTime(tsa):tsa.dt:EndTime(tsa);
R = R';
if nargin == 2
   switch (tsflag)
   case 'sec'
      R = R;
   case 'sec0'
      R = (R - StartTime(tsa));     
   case 'ts'
	  R = R*10000;
   case 'ms'
      R = R*1000;
   otherwise
      error('Range called with invalid tsflag.');
   end
end


