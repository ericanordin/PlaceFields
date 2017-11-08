function [rotatedMatrix, nan_map ] = rotateAndPrep(unrotatedMatrix, angle, trimLength, cornerCompare)
%rotateAndPrep rotates and prepares a matrix for comparison
%   imrotate fills in new spaces with '0'. This function converts those
%   additions to NaNs, trims the rotatedMatrix down to its original size so
%   that it can be compared to the standard matrix, marks the location of
%   the NaN values, and converts all NaN values to 0 in preparation of
%   matrix correlation.

% Input:
% unrotatedMatrix = the matrix before rotation
% angle = CCW angle in degrees that unrotatedMatrix will be rotated by
% trimLength = the length of unrotatedMatrix which rotatedMatrix will be
% trimmed to
% cornerCompare = the sample matrix which marks the location of the corners
% for the matrix following rotation

% Output
% rotatedMatrix: the matrix following rotation and NaN conversion
% nan_map: the matrix which marks the locations of the NaN values in
% rotatedMatrix

%Rotate
rotatedMatrix = imrotate(unrotatedMatrix, angle, 'bilinear');


%Convert regions outside of buffer to NaN

%Determine boundaries of cut_map portion of
%buffered_map
[cornerRow, cornerColumn] = find(cornerCompare);
leftLimit = min(cornerColumn)-1;
rightLimit = max(cornerColumn)+1;
lowerLimit = max(cornerRow)+1;
upperLimit = min(cornerRow)-1;

%Fill in 0s added by imrotate with NaNs by making
%everything past boundaries of cut_map NaN
rotatedMatrix(1:leftLimit,:) = NaN;
rotatedMatrix(rightLimit:end,:) = NaN;
rotatedMatrix(:,1:upperLimit) = NaN;
rotatedMatrix(:,lowerLimit:end) = NaN;

%The amount to be trimmed off of each side
%Assumes square trimming
trimSide = (size(rotatedMatrix,1) - trimLength)/2;
%Must account for whole (remainer == 0) or non-whole
%(remainder ~= 0) trimSide value for even or odd size
%of temp_map
remainder = rem(trimSide, floor(trimSide));
if remainder == 0
    %Same trim taken from all sides
    leftBound = trimSide;
    upperBound = trimSide;
    rightBound = size(rotatedMatrix,2) - trimSide;
    lowerBound = size(rotatedMatrix,1) - trimSide;
else
    if leftLimit > (size(rotatedMatrix,2) - rightLimit)
        %More trim taken from left buffer
        leftBound = ceil(trimSide);
        rightBound = size(rotatedMatrix,2) - floor(trimSide);
    else
        %More trim taken from right buffer
        leftBound = floor(trimSide);
        rightBound = size(rotatedMatrix,2) - ceil(trimSide);
    end
    
    if upperLimit > (size(rotatedMatrix,1) - lowerLimit)
        %More trim taken from upper buffer
        upperBound = ceil(trimSide);
        lowerBound = size(rotatedMatrix,1) - floor(trimSide);
    else
        %More trim taken from lower buffer
        upperBound = floor(trimSide);
        lowerBound = size(rotatedMatrix,1) - ceil(trimSide);
    end
end

rotatedMatrix = rotatedMatrix(upperBound+1:lowerBound, ...
    leftBound+1:rightBound);
nan_map = isnan(rotatedMatrix(:,:));
rotatedMatrix(nan_map) = 0;  % turn all NaN's to zero


end

