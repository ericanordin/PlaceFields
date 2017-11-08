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



end

