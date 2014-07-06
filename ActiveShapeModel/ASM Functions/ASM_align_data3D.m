function [Vertices,tform]=ASM_align_data3D(Vertices,VerticesB)
% Aligns the contours positions, center the data and remove rotation

% Center data to remove translation 
offsetVertices = -mean(Vertices,1);
Vertices = Vertices + repmat(offsetVertices,size(Vertices,1),1);

offsetVerticesB = -mean(VerticesB,1);
VerticesB = VerticesB + repmat(offsetVerticesB,size(VerticesB,1),1);

% Correct for rotation
% Calculate angle to center of all points
rot = atan2(Vertices(:,2),Vertices(:,1));
rotB = atan2(VerticesB(:,2),VerticesB(:,1));

% Subtract the mean angle
offsetrxy=mean(rotB-rot);
rot = rot + offsetrxy;

% Make the new points, which all have the same rotation
dist = sqrt(Vertices(:,1).^2+Vertices(:,2).^2);
Vertices(:,1) =dist.*cos(rot);
Vertices(:,2) =dist.*sin(rot);

% Correct for rotation
% Calculate angle to center of all points
rot = atan2(Vertices(:,3),Vertices(:,2));
rotB = atan2(VerticesB(:,3),VerticesB(:,2));

% Subtract the mean angle
offsetryz=mean(rotB-rot);
rot = rot + offsetryz;

% Make the new points, which all have the same rotation
dist = sqrt(Vertices(:,2).^2+Vertices(:,3).^2);
Vertices(:,2) =dist.*cos(rot);
Vertices(:,3) =dist.*sin(rot);

% Store transformation object
tform.offsetV=offsetVertices;
tform.offsetrxy=offsetrxy;
tform.offsetryz=offsetryz;

