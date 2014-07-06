function V=ASM_align_data_inverse3D(V,tform)

rot = atan2(V(:,3),V(:,2));
rot = rot - tform.offsetryz;

dist = sqrt(V(:,2).^2+V(:,3).^2);
V(:,2) =dist.*cos(rot);
V(:,3) =dist.*sin(rot);

rot = atan2(V(:,2),V(:,1));
rot = rot - tform.offsetrxy;

dist = sqrt(V(:,1).^2+V(:,2).^2);
V(:,1) =dist.*cos(rot);
V(:,2) =dist.*sin(rot);

V = V - repmat(tform.offsetV,size(V,1),1);

