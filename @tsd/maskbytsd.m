function C = MaskByTSD(A,B)
%
% C = ts/MaskByTSD(A,B)
% 
% INPUTS:
%    A = tsd object
%    B = tsd object
%
% OUTPUTS:
%    C = tsd only including times in which nearest tsd sample from B
%          is not NaN
%
% ADR 2000
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0.
% Version control M3.0.

TIMES = Range(A, 'ts');
DATAS = Data(A);
nT = length(TIMES);

tTime = Range(B, 'ts');
tData = Data(B);

for iS = 1:nT
  if isnan(tData(binsearch(tTime, TIMES(iS))))
    TIMES(iS) = nan;
  end
end

f = find(~isnan(TIMES));
C = tsd(TIMES(f), DATAS(f));
