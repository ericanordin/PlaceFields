function T1 = EndTime(tsa, tsflag)
%
% T1 = tsd/EndTime(tsd, tsflag)
%	returns last timestamp covered by tsa
%
%      tsflag: if 'ts' returns time in timestamps (default),
%              if 'sec' returns time in sec
%              if 'ms' returns time in ms
%
%
%
% ADR 1998
% version L4.0
% status PROMOTED
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0.
% Version control M3.0.

T1 = max(tsa.t);

if nargin == 2
   switch tsflag
   case 'sec'
      T1 = T1/10000;
   case 'ms'
      T1 = T1/10;
   case 'ts'
   otherwise
      error('Unknown tsflag.');
   end
end
