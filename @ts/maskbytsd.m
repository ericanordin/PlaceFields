function ts1 = MaskByTSD(ts0, tsd0)
%
% ts1 = ts/MaskByTSD(ts0, tsd0)
% 
% INPUTS:
%    ts0 = ts object
%    tsd0 = tsd object
%
% OUTPUTS:
%    ts1 = ts only including times in which nearest tsd sample
%          is not NaN
%
% ADR 2000
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0.
% Version control M3.0.

SpikesIn = Data(ts0);
nSpikes = length(SpikesIn);

tTime = Range(tsd0, 'ts');
tData = Data(tsd0);

for iS = 1:nSpikes
  if isnan(tData(binsearch(tTime, SpikesIn(iS))))
    SpikesIn(iS) = nan;
  end
end

ts1 = ts(SpikesIn(find(~isnan(SpikesIn))));
