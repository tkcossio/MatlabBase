function [gtc,dgtc]=ASM_getProfileAndDerivatives2D(I,P,N,k)
% This function getProfileAndDerivatives, samples the image on lines
% perpendicular to the contour points, which are used to train the
% appearance model, and to find the best location for a contour point later
% on.
%
% [gt,dgt]=getProfileAndDerivatives(I,P,N,k)
%
% inputs,
%  I : The image color or greyscale
%  P : The locations of the contour points
%  N : The normals of the contour points
%  k : The length of the sampled line in one normal direction, 
%      total length becomes k*2+1
%
% outputs,
%  gt : Columns with the sampled lines perpendicular to the contour
%  dgt : The first order derivatives of the sampled lines
%
% Function written by D.Kroon University of Twente (February 2010)



gtc =zeros((k*2+1)*size(I,3),size(P,1));
dgtc=zeros((k*2+1)*size(I,3),size(P,1));

for i=1:size(I,3)
    xi=linspace_multi(P(:,1)-N(:,1)*k,P(:,1)+N(:,1)*k,k*2+1);
    yi=linspace_multi(P(:,2)-N(:,2)*k,P(:,2)+N(:,2)*k,k*2+1);
    xi(xi<1)=1; xi(xi>size(I,1))=size(I,1);
	yi(yi<1)=1; yi(yi>size(I,2))=size(I,2);

    % Sample on the normal lines
    gt= interp2(I(:,:,i)',xi,yi,'cubic')';
    
    gt(isnan(gt))=0;
    % Get the derivatives
    dgt=[gt(2,:)-gt(1,:);(gt(3:end,:)-gt(1:end-2,:))/2;gt(end,:)-gt(end-1,:)];
    % Store the grey profiles and derivatives for the different color
    % channels
    b=(k*2+1)*(i-1)+1; e=b+(k*2+1)-1;
    gtc(b:e,:)=gt; dgtc(b:e,:)=dgt;
end
% Normalize the derivatives
dgtc=dgtc./repmat((sum(abs(dgtc),1)+eps),(k*2+1)*size(I,3),1);
%gtc=gtc./repmat((sum(abs(gtc),1)+eps),(k*2+1)*size(I,3),1);