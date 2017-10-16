function mTS = Mask(ts, varargin)
%
% mTS = ts/Mask(ts, TrialPairs....)
%
% INPUTS:
%    ts = ts object
%    TrialPairs = pairs of start/end times (can be matrices of n x 2)
%
% OUTPUTS:
%    mts = ts object with times *not* in TrialPairs set to NaN
%
% SLC 1999
% Modified the tsd Mask method to work for ts objects.
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0.
% Version control M3.0.

StatusWarning('UNKOWN', '@ts/Mask');

% Unwrap trial pairs
MaskOFF = [StartTime(ts)-1];
MaskON = []; 
for iTP = 1:length(varargin)
   curMask = varargin{iTP};
   MaskOFF = cat(1, MaskOFF, curMask(:,2));
   MaskON = cat(1,MaskON, curMask(:,1));
end
MaskON = cat(1, MaskON, EndTime(ts)+1);
nTransitions = length(MaskON);
MaskON = sort(MaskON);
MaskOFF = sort(MaskOFF);

% Construction output ts
mTS = ts;

% Now implement mask
for iT = 1:nTransitions
   f = find(mTS.t > MaskOFF(iT) & mTS.t < MaskON(iT));
   mTS.t(f) = NaN;
   %mTSD.data(f) = NaN;
end
