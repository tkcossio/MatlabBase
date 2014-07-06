function [Vertices,tform]=ASM_align_data2D(Vertices)
% Aligns the contours positions, center the data and remove rotation

% Center data to remove translation 
offsetv = -mean(Vertices,1);
Vertices(:,1) = Vertices(:,1) + offsetv(1);
Vertices(:,2) = Vertices(:,2) + offsetv(2);

%offsets=mean(sqrt(Vertices(:,1).^2+Vertices(:,2) .^2));
%Vertices(:,1)=Vertices(:,1)/offsets;
%Vertices(:,2)=Vertices(:,2)/offsets;

% Correct for rotation
% Calculate angle to center of all points
rot = atan2(Vertices(:,2),Vertices(:,1));
% Subtract the mean angle
offsetr=-mean(rot(1:round(end/2)));
rot = rot+offsetr;
% Make the new points, which all have the same rotation
dist = sqrt(Vertices(:,1).^2+Vertices(:,2).^2);
Vertices(:,1) = dist.*cos(rot);
Vertices(:,2) = dist.*sin(rot);

% Store transformation object
tform.offsetv=offsetv;
tform.offsetr=offsetr;
%tform.offsets=offsets;

