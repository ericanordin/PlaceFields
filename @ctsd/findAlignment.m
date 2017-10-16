function ix = findAlignment(D, tstmp)

% tsArray/findAlignment
% 	ix = findAlignment(D, tstmp)
%
% 	Returns an index i= such that D(ix) occurs 
%		at timestamp tstmp

% ADR
% version: L4.0
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0.
% Version control M3.0.

ix = round((tstmp - D.t0)/D.dt) + 1;