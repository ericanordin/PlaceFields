function ts = findTime(D, ix)
%
% ctsd/findTime
% 	ix = findTime(D, ix)
%
% 	Returns timestamsp at which index ix occues

% ADR
% version L4.0
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0.
% Version control M3.0.

ts = D.t0 + ix * D.dt;
