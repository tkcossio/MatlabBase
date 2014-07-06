function Vertices=AAM_align_data_inverse2D(Vertices,tform)

tform.offsets=sqrt( tform.offsetsx^2+tform.offsetsy^2);
tform.offsetr=atan2(tform.offsetsy  ,tform.offsetsx);

% Correct for rotation
rot = atan2(Vertices(:,2),Vertices(:,1));
rot = rot-tform.offsetr;
dist = sqrt(Vertices(:,1).^2+Vertices(:,2).^2);

Vertices(:,1) =dist.*cos(rot);
Vertices(:,2) =dist.*sin(rot);

Vertices = Vertices / tform.offsets;
Vertices(:,1) = Vertices(:,1) - tform.offsetv(1);
Vertices(:,2) = Vertices(:,2) - tform.offsetv(2);




