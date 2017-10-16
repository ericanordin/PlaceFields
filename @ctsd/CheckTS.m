function bool = CheckTS(varargin)

% bool = [c]tsd/CheckTS(X0, X1, X2, ...)
%
% checks to make sure that all timestamps are identical for all tsds included.
% works with combinations of ctsd and tsd

% ADR 1998
% version L4.0
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0.
% Version control M3.0.

adrlib;

% for speed, first check to make sure if they're all ctsd's...
flagAllCTSD = true;
for iX = 1:length(varargin)
   if ~isa(varargin{iX}, 'ctsd')
      flagAllCTSD = false;
   end
end

if flagAllCTSD
   % if they're all ctsd's, just check DT's and start/end times
   DT0 = DT(varargin{1});
   T0 = StartTime(varargin{1});
   T1 = EndTime(varargin{1});
   for iX = 2:length(varargin)
      if (DT0 ~= DT(varargin{iX}))
         bool = false; 
         return
      elseif (T0 ~= StartTime(varargin{iX}))
         bool = false;
         return
      elseif (T1 ~= EndTime(varargin{iX}))
         bool = false;
         return
      end % if
   end % for
   bool = true;
   return
else
   R0 = Range(varargin{1});
   for iX = 2:length(varargin)
      R1 = Range(varargin{iX});
      if (length(R0) ~= length(R1))
         bool = false;                  % if not same length, can't be equal.
         return
      end
      if (min(R0 == R1) == 0)
         bool = false;                  % if there are any non-equal elts, not equal
         return
      end
   end
   bool = true;                        % nothing failed, must be ok
   return
end
