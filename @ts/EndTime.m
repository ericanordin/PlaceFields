function t = EndTime(TS, tsflag)

% t = ts/EndTime(TS. tsflag)
%
% returns last timestamp in TS
%      tsflag: if 'ts' returns time in timestamps (default),
%              if 'sec' returns time in sec
%              if 'ms' returns time in ms
%
%
% ADR
% version L4.0
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0.
% Version control M3.0.

t = max(TS.t);

if nargin == 2
   switch tsflag
   case 'sec'
      t = t/10000;
   case 'ms'
      t = t/10;
   case 'ts'
   otherwise
      error('Unknown tsflag.');
   end
end
