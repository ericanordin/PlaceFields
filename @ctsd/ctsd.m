function tsa = ctsd(qT0, qDT, qData)
%
% tsa = ctsd(qT0, qDT, qData)
%
% ctsd is a class of "timestamped arrays"
% 	It includes a starting time (T0), a deltaT (DT)
% 	and data (possibly an array).  Data begins at T0
%	and proceeds sequentially at a timestep of DT.
% 	It is completely compatible with tsd.
%    If Data is n-dimensional, then time should be the FIRST axis.
%
% Methods
%    ctsd/Range     - Timestamps used
%    ctsd/Data      - Returns the data component
%    ctsd/DT        - Returns the DT value
%    ctsd/StartTime - First timestamp
%    ctsd/EndTime   - Last timestamp
%    ctsd/Restrict  - Keep data within a certain range
%    ctsd/CheckTS   - Makes sure that a set of tsd & ctsd objects have identical start and end times
%    ctsd/cat       - Concatenate ctsd and tsd objects
%    ctsd/Mask      - Make all non-mask values NaN

% ADR 
% version L4.0
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0.
% Version control M3.0.


tsa.t0 = NaN;
tsa.dt = NaN;
tsa.data = NaN;

switch nargin
   
case 0
   tsa.t0 = NaN;
   tsa.dt = NaN;
   tsa.data = NaN;
case 1
   if isa(qT0, 'tsd')
      % convert tsd to ctsd
      tsa.t0 = StartTime(qT0);
      tsa.dt = median(diff(Range(qT0,'ts')));
      tsa.data = interp1q(Range(qT0,'ts'), Data(qT0), (StartTime(qT0):tsa.dt:EndTime(qT0))');
   elseif isa(qT0, 'ctsd')
      tsa = qT0;
      return;
   else
      error('Unknown copy-from object');
   end
case 3
   tsa.t0 = qT0;
   tsa.dt = qDT;
   tsa.data = qData;
otherwise
   error('Constructor error: ctsd');
end
tsa = class(tsa, 'ctsd');
