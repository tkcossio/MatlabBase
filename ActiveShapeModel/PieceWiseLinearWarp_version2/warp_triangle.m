function J=warp_triangle(I,xy,uv,ImageSize,tri)
% This function WARP_TRIANGLE will perform a piecewise linear 2D image warp 
% base ond control point pairs. The function is almost equal to cp2tform 
% 'piecewise linear' in combination with imtransform 'bicubic'. The pixel 
% interpolation used is bi-cubic interpolation, and boundary condition clamp.
%
%   J=warp_triangle(I,xy,uv,ImageSize,tri)
%
% inputs,
%   I : 2D color or greyscale input image
%   xy : M x 2, Control point coordinates in the input image
%   uv : M x 2, Control point coordinates in the warped/output image
% (optional)
%   ImageSize : 1 x 2 sizes of warped/output image, or set to initial output 
%               image.
%   tri : N x 3 Control point triangulation list, suchs as produced by delaunay
%
% outputs,
%   J : The warped image
%
%
% note:
%  C-coded files are available for speed. Run compile_c_files to build the
%  fast Mex files
%
% example,
%  I = im2double(imread('images/lena.bmp'));
%  xy=[1 1;1 size(I,2);size(I,1) 1;size(I,1) size(I,2)];
%  uv=[256 512; 1 256; 512  256; 256 1];
%
%  ImageSize=[512 512];
%  J = warp_triangle(I,xy,uv,ImageSize);
%  figure,
%  subplot(1,2,1), imshow(I)
%  subplot(1,2,2), imshow(J)
%
%    
% Function is written by D.Kroon University of Twente (July 2011)

functionname='warp_triangle.m';
functiondir=which(functionname);
functiondir=functiondir(1:end-length(functionname));
addpath([functiondir '/functions'])

if(nargin<4), ImageSize=[size(I,1) size(I,2)]; end
if(nargin<5), tri = delaunay(xy(:,1),xy(:,2)); end

if((size(xy,2)~=2)||(size(uv,2)~=2))
    error('warp_triangle:inputs', 'uv and xy must have sizes m x 2');
end
if(size(xy,1)~=size(uv,1))
    error('warp_triangle:inputs', 'uv and xy must have the same size');
end

% Matlab convention
xy=xy(:,[2 1]);
uv=uv(:,[2 1]);

classI=class(I);
if(~strcmpi(classI,'double'))
    I=double(I);
    if(numel(ImageSize)>3), ImageSize=double(ImageSize); end
end
if(numel(ImageSize)>3)
    J=warp_triangle_double(I,double(xy),double(uv),double(tri),double(size(ImageSize)),ImageSize);
else
    J=warp_triangle_double(I,double(xy),double(uv),double(tri),double(ImageSize));
end

if(~(strcmpi(classI,'double')||strcmpi(classI,'single')))
    J=cast(J,classI);
end

