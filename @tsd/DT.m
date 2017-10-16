function dt = DT(tsa, tsflag)
%
% dt = tsd/DT(tsd, tsflag)
%	returns DT from tsa
% INPUTS
%      tsd 
%      tsflag: if 'ts' returns time in timestamps (default),
%              if 'sec' returns time in sec
%              if 'ms' returns time in ms
%
% ADR 
% version L1.0
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0.
% Version control M3.0.

dt = mean(diff(tsa.t));

if nargin == 2
   switch tsflag
   case 'sec'
      dt = dt/10000;
   case 'ms'
      dt = dt/10;
   case 'ts'
   otherwise
      error('Unknown tsflag.');
   end
end
