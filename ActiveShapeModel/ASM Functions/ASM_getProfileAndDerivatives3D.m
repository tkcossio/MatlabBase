function [gt,dgt]=ASM_getProfileAndDerivatives3D(I,pV,nV,k)
% This function getProfileAndDerivatives, samples the image on lines
% perpendicular to the contour points, which are used to train the
% appearance model, and to find the best location for a contour point later
% on.
%
% [gt,dgt]=getProfileAndDerivatives(I,pV(:,1),pV(:,2),nV(:,1),nV(:,2),k)
%
% inputs,
%  I : The image color or greyscale
%  pV: The locations of the contour points
%  nV: The normals of the contour points
%  k : The length of the sampled line in one normal direction,
%      total length becomes k*2+1
%
% outputs,
%  gt : Columns with the sampled lines perpendicular to the contour
%  dgt : The first order derivatives of the sampled lines
%
% Function written by D.Kroon University of Twente (February 2010)

xi=linspace_multi(pV(:,1)-nV(:,1)*k,pV(:,1)+nV(:,1)*k,k*2+1);
yi=linspace_multi(pV(:,2)-nV(:,2)*k,pV(:,2)+nV(:,2)*k,k*2+1);
zi=linspace_multi(pV(:,3)-nV(:,3)*k,pV(:,3)+nV(:,3)*k,k*2+1);

xi(xi<1)=1; xi(xi>size(I,1))=size(I,1);
yi(yi<1)=1; yi(yi>size(I,2))=size(I,2);
zi(zi<1)=1; zi(zi>size(I,3))=size(I,3);


% Sample on the normal lines
%gt= interp3(permute(I,[2 1 3]),xi,yi,zi,'cubic')';
gt= interpfast(I,xi,yi,zi,'cubic')';

%gt(isnan(gt))=0;

% Get the derivatives
dgt=[gt(2,:)-gt(1,:);(gt(3:end,:)-gt(1:end-2,:))/2;gt(end,:)-gt(end-1,:)];

% Store the grey profiles and derivatives for the different color
% channels

% Normalize the derivatives
if(nargout>1)
    dgt=dgt./repmat((sum(abs(dgt),1)+eps),k*2+1,1);
end
%gt=gt-mean(gt(:));
%gt=gt./std(gt(:));
