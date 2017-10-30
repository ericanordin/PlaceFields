%From file exchange. Will be modified.

% This function takes an image or matrix and rotate it by theta (radian, rigid motion)
% The image could also be enlarged or shrinked by an optional scaling factor
%
% img:           input image (or or a 2D Matrix)
% dtheta:        the amount of rotation in radian
% scale:         a scalar number that scales the image (default=1)
% interpMethod:  
%      'nearest' - nearest neighbor interpolation
%      'linear'  - bilinear interpolation (defalt)
%      'spline'  - spline interpolation
%      'cubic'   - bicubic interpolation
% padding:       background filling color (0~255, black~white). Default: black
%
% Example:
%   im = imread('cameraman.tif')
%   imshow(imTransform(im,pi/6)) 

function imT = imTransform(img,dtheta,scale,interpMethod,padding)
if nargin<5
    padding=0;
end

if nargin<4
    interpMethod='linear'; %Turn of repeated edge rendering
end

if nargin<3
    scale=1;
end

if ~isfloat(img)
    u8=true;
    img=double(img); %Make sure the operation is done for floats
end
ND=ndims(img);
if ND==3,
    [yDim,xDim,~] = size(img); %Obtain the number of rows and columns
else
    [yDim,xDim] = size(img);
end

rotMat=[cos(dtheta),sin(dtheta);-sin(dtheta),cos(dtheta)];
xMax=(xDim-1)*scale;
yMax=(yDim-1)*scale;
[xx0,yy0]=meshgrid(0:xMax,0:yMax);

ctr=[xMax/2;yMax/2]; %Find the center location for rotation operation

% Build a pseudo matt for the target image
corners=[0,0,xMax,xMax;0,yMax,0,yMax];
corners(1,:)=corners(1,:)-ctr(1);
corners(2,:)=corners(2,:)-ctr(2);
cornersT=rotMat*corners;
cornerMin=min(cornersT,[],2);
cornerMax=max(cornersT,[],2)-cornerMin+1;%Maximum pixel after transform
xyNew=ceil(cornerMax)-1;
[xx,yy]=meshgrid(0:xyNew(1),0:xyNew(2));

%Map back to the original image
xx=xx+cornerMin(1);
yy=yy+cornerMin(2);

%Rotate backwards and shift to (original location -1),Finally scale down and add 1
xxMap = (rotMat(1)*xx-rotMat(3)*yy+ctr(1))/scale;
yyMap = (-rotMat(2)*xx+rotMat(4)*yy+ctr(2))/scale;


if ND==3
    imT(:,:,1)=interp2(img(:,:,1), xxMap,yyMap,interpMethod);
    imT(:,:,2)=interp2(img(:,:,2), xxMap,yyMap,interpMethod);
    imT(:,:,3)=interp2(img(:,:,3), xxMap,yyMap,interpMethod);
else
    imT=interp2(img, xxMap,yyMap,interpMethod);
end

imT(isnan(imT)) = padding; %Fill black pixels
% 
if u8
    imT=uint8(imT);
end

end
