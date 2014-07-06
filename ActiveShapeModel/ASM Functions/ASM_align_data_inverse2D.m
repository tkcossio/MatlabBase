function pos=ASM_align_data_inverse2D(pos,tform)
% Correct for rotation
rot = atan2(pos(:,2),pos(:,1));
rot = rot-tform.offsetr;
dist = sqrt(pos(:,1).^2+pos(:,2).^2);
% *tform.offsets;
pos(:,1) =dist.*cos(rot);
pos(:,2) =dist.*sin(rot);
pos(:,1) = pos(:,1) - tform.offsetv(1);
pos(:,2) = pos(:,2) - tform.offsetv(2);
