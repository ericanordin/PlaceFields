function tsa = ts(t)
%
% TSD = ts(t)
%
% A ts object contains a sequence of timestamps as
% might be stored in an NSMA t-file.
%
% Methods
%    ts/Data         - Returns the timestamps as a matlab array
%    ts/StartTime    - First timestamp
%    ts/EndTime      - Last timestamp
%    ts/Restrict      - Keep timestamps within a certain range
%
% ADR
% version L4.1
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0.
% Version control M3.0.

if nargin == 0
   tsa.t = [];
   return
end

tsa.t = t;
tsa = class(tsa, 'ts');
