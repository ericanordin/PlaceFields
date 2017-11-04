%temp_map has 0s substituted for NaNs after NaNs were counted in cut_map.
%tempbuf adds a buffer of NaN values around temp_map, so after rotation all
%0s surrounding the buffer can be converted to NaN.
%May make a larger buffer to prevent data loss during rotation.
tempbuf = NaN(size(temp_map,1)+2, size(temp_map,2)+2);
tempbuf(2:end-1, 2:end-1) = temp_map;
t2 = imrotate(tempbuf, 30);
imagesc(t2);
%nan_i may need buffer and rotation as well

%Buffer size should be diagonal length of image so that rotation will never
%go out of bounds
%51x51 = 72.12 diagonal -> 73
%Trim down image to be size of unrotated image after rotation has occurred
