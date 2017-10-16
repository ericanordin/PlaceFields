function R = Restrict(D, t0, t1)

% ctsd/Restrict
% 	R = Restrict(tsa, t0, t1)
% 	Returns a new tsa (ctsd) R so that D.Data is between 
%		timestamps t0 and t1
%
%    R = Restrict(tsa, t)
%    Returns a new tsd (not ctsd) R so that D.Data includes
%         only those timestamps in t

% ADR 
% version L5.0
% v4.1 29 oct 1998 now can handle nargin=2
% v5.0 30 oct 1998 time dimension is always 1st dimension
% v5.1 19 jan 1998 now can handle t0 and t1 as arrays
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0.
% Version control M3.0.

switch nargin
case 2                             % R = Restrict(tsd, t)
   ix = findAlignment(D, t0);
   tlist = Range(D, 'ts');
   R = tsd(tlist(ix), SelectAlongFirstDimension(D.data, ix));
   
case 3                             % R = Restrict(tsd, t0, t1)
   if length(t0) ~= length(t1)
      error('t0 and t1 must be the same length.');
   end
   if length(t0) == 1
      R0 = max(1, findAlignment(D, t0));
      R1 = min(length(D.data), findAlignment(D, t1));
      
      R = ctsd(t0, D.dt, SelectAlongFirstDimension(D.data, R0:R1));
   else
      D.data = full(D.data);
      DATA = [];
      TIME = [];
      for it = 1:length(t0)
         R0 = max(1, findAlignment(D, t0(it)));
         R1 = min(length(D.data), findAlignment(D, t1(it))); 
         TIME = cat(2, TIME, findTime(D, R0):D.dt:findTime(D, R1));
         DATA = cat(1, DATA, SelectAlongFirstDimension(D.data, R0:R1));
      end
      R = tsd(TIME', DATA);
   end
   
otherwise
   error('Unknown number of input arguments.');
   
end % switch
