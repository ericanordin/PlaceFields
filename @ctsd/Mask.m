function mTSD = Mask(tsd, varargin)
%
% mTSD = ctsd/Mask(tsd, TrialPairs....)
%
% INPUTS:
%    tsd = tsd object
%    TrialPairs = pairs of start/end times (can be matrices of n x 2)
%
% OUTPUTS:
%    mTSD = ctsd object with times *not* in TrialPairs set to NaN

% ADR 1998
% version L4.0
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0.
% Version control M3.0.

% Unwrap trial pairs
MaskOFF = [StartTime(tsd)-1];
MaskON = []; 
for iTP = 1:length(varargin)
   curMask = varargin{iTP};
   MaskOFF = cat(1, MaskOFF, curMask(:,2));
   MaskON = cat(1,MaskON, curMask(:,1));
end
MaskON = cat(1, MaskON, EndTime(tsd)+1);
nTransitions = length(MaskON);
MaskON = sort(MaskON);
MaskOFF = sort(MaskOFF);

% Construction output tsd
mTSD = tsd;

% Now implement mask
for iT = 1:nTransitions
   m0 = findAlignment(tsd, MaskOFF(iT));
   m1 = min(length(tsd.data), findAlignment(tsd, MaskON(iT)));
   if m0 <= 1; 
      m0 = 1;
   elseif findTime(tsd, m0) == MaskOFF(iT)
      m0 = m0 + 1;
   end       % really want > not >=
   if m1 >= length(tsd.data)
      m1 = length(tsd.data);
   elseif findTime(tsd, m1) == MaskON(iT)
      m1 = m1 - 1;
   end       % really want < not <=
   mTSD.data(m0:m1) = NaN;
end
