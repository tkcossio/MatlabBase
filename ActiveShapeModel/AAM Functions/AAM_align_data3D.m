function [Vertices,tform]=AAM_align_data3D(Vertices,VerticesB,options)
% Remove rotation and translation and scale : Procrustes analysis 


% Center data to remove translation 
offsetv = -mean(Vertices,1);
Vertices(:,1) = Vertices(:,1) + offsetv(1);
Vertices(:,2) = Vertices(:,2) + offsetv(2);
Vertices(:,3) = Vertices(:,3) + offsetv(3);

offsetvB = -mean(VerticesB,1);
VerticesB(:,1) = VerticesB(:,1) + offsetvB(1);
VerticesB(:,2) = VerticesB(:,2) + offsetvB(2);
VerticesB(:,3) = VerticesB(:,3) + offsetvB(3);

% Warp-Matrix between point sets
S=Vertices'*VerticesB;
M=[(S(1,1)+S(2,2)+S(3,3))  (S(2,3)-S(3,2))      (S(3,1)-S(1,3))      (S(1,2)-S(2,1));...
   (S(2,3)-S(3,2))      (S(1,1)-S(2,2)-S(3,3))  (S(1,2)+S(2,1))      (S(3,1)+S(1,3));...
   (S(3,1)-S(1,3))      (S(1,2)+S(2,1))     (-S(1,1)+S(2,2)-S(3,3))  (S(2,3)+S(3,2));...
   (S(1,2)-S(2,1))      (S(3,1)+S(1,3))      (S(2,3)+S(3,2))      (-S(1,1)-S(2,2)+S(3,3))];

%Compute eigenvalues
[E,~] = eig(M);

% Calculate the quaternion
quat = E(:,4);
[~,ind]=max(abs(quat)); 
quat=quat.*sign(quat(ind)); 
quat=quat./norm(quat);
offsetq=quat;

%Compute the rotation matrix
q0=quat(1); qx=quat(2);  qy=quat(3); qz=quat(4); qxyz =quat(2:4);
R=qxyz*qxyz.' + [q0 -qz qy; qz q0 -qx; -qy qx  q0]^2;

% Calculate Size
VerticesR=(R*Vertices')';
if(options.scale3)
    offsets(1)=sum(VerticesB(:,1).*VerticesR(:,1))/sum(VerticesR(:,1).^2);
    offsets(2)=sum(VerticesB(:,2).*VerticesR(:,2))/sum(VerticesR(:,2).^2);
    offsets(3)=sum(VerticesB(:,3).*VerticesR(:,3))/sum(VerticesR(:,3).^2);
    VerticesR(:,1)=VerticesR(:,1)*offsets(1);
    VerticesR(:,2)=VerticesR(:,2)*offsets(2);
    VerticesR(:,3)=VerticesR(:,3)*offsets(3);
else
    summ = @(M) sum(M(:));
    offsets=summ( VerticesB.*VerticesR)/summ(VerticesR.^2);
    VerticesR=VerticesR*offsets;
end
Vertices=VerticesR;
   
% Store transformation object
tform.offsetv=offsetv;
if(options.scale3)
    tform.offsetq=offsetq; 
else
    tform.offsetq=offsetq*offsets;
end
tform.offsets=offsets;
tform.offsetq=tform.offsetq*options.scaleqmatrix;

