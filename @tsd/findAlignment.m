function ix = findAlignment(D, tstmp)

% tsArray/findAlignment
% 	ix = findAlignment(D, tstmp)
%
% 	Returns an index i= such that D(ix) occurs 
%		at timestamp tstmp
% 	Finds closest index.
%
% ADR
% version: L4.0
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0.
% Version control M3.0.

% modified 13 Sep 02 by ncst to use nearest neighbor interpolation

nIX = length(tstmp);

% ix = zeros(size(tstmp));
% for iIX = 1:nIX
%    ix(iIX) = binsearch(D.t, tstmp(iIX));
% end

ix = interp1(D.t,1:length(D.t),tstmp,'nearest');
% out of range values (tstmp values outside the range of D.t) are returned as NaNs, 
% so remove the NaNs before returning
    ix = ix(find(~isnan(ix)));
