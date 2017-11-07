[X,Y] = meshgrid(-51:51, -51:51);
%Produces indeces relative to 0 (range from -51 to 51)
circle_i = sqrt(X.^2 + Y.^2)<5; %5 pixels on either side of the centre point
%Index of 0s and 1s

imagesc(circle_i); %Creates circle
temp_map(circle_i) = 0; %Create zeros in central circle where circle_i == 1
%Only values corresponding to circle_i == 1 are impacted