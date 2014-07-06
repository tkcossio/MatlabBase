function J=warp_tetrahedron(I,xyz,uvw,ImageSize,tetra)
% This function WARP_TETRAHEDRON will perform a piecewise 3D linear image 
% warp base ond control point pairs. The pixel interpolation used is 
% tri-linear interpolation, and boundary condition clamp.
%
%   J=warp_tetrahedron(I,xyz,uvw,ImageSize,tetra)
%
% inputs,
%   I : 3D greyscale input image
%   xyz : M x 3, Control point coordinates in the input image
%           (Matlab Y,X,Z convention is used)
%   uvw : M x 3, Control point coordinates in the warped/output image
% (optional)
%   ImageSize : 1 x 3 sizes of warped/output image, or set to initial output 
%               image volume
%   tetra : N x 4 Control point triangulation list, suchs as produced by
%   delaunayn
%
% outputs,
%   J : The warped 3D image
%
% note:
%  C-coded files are available for speed. Run compile_c_files to build the
%  fast Mex files
%
% example,
%  load('images/testdata.mat');
%  I=single(I);
%  uvw=[1 1                1;
%      1 size(I,2)         1;
%      size(I,1) 1         1;
%      size(I,1) size(I,2) 1;
%      1 1                 size(I,3);
%      1 size(I,2)         size(I,3);
%      size(I,1) 1         size(I,3);
%      size(I,1) size(I,2) size(I,3)];
%  xyz=[(uvw(:,2)+uvw(:,1))/1.5-32 (uvw(:,2)-uvw(:,1))/1.5+32 uvw(:,3)]; 
%  ImageSize=[64 64 64];
%
%  J=warp_tetrahedron(I,xyz,uvw,ImageSize);
%  figure,
%  subplot(2,3,1), imshow(squeeze(I(32,:,:)));
%  subplot(2,3,2), imshow(squeeze(I(:,32,:)));
%  subplot(2,3,3), imshow(I(:,:,32));
%  subplot(2,3,4), imshow(squeeze(J(32,:,:)));
%  subplot(2,3,5), imshow(squeeze(J(:,32,:)));
%  subplot(2,3,6), imshow(J(:,:,32));
%    
% Function is written by D.Kroon University of Twente (July 2011)

xyz=xyz(:,[2 1 3]);
uvw=uvw(:,[2 1 3]);

if(nargin<4), 
    ImageSize=[size(I,1) size(I,2) size(I,3)];
end
if(nargin<5),
    tetra= delaunayn(xyz,{'Qt','Qbb','Qc','Qz'});
end

classI=class(I);
if(~(strcmpi(classI,'double')||strcmpi(classI,'single')))
    I=single(I);
    if(numel(ImageSize)>3), ImageSize=single(ImageSize); end
end
if(numel(ImageSize)>3)
    J=warp_tetrahedron_double(I,double(xyz),double(uvw),double(tetra),size(ImageSize),ImageSize);
else
    J=warp_tetrahedron_double(I,double(xyz),double(uvw),double(tetra),double(ImageSize));
end

if(~(strcmpi(classI,'double')||strcmpi(classI,'single')))
    J=cast(J,classI);
end



 