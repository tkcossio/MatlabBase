function Vertices=AAM_align_data_inverse3D(Vertices,tform,options)
tform.offsetq=tform.offsetq/options.scaleqmatrix;
qn=norm(tform.offsetq,2);
quat=tform.offsetq/qn;
if(options.scale3)
    offsets=tform.offsets;
else
    offsets=qn;
end
offsetv=tform.offsetv;

%Compute the rotation matrix
q0=quat(1); qx=quat(2);  qy=quat(3); qz=quat(4); qxyz =quat(2:4);
R=qxyz*qxyz.' + [q0 -qz qy; qz q0 -qx; -qy qx  q0]^2;
invR=inv(R);

if(length(offsets)>1)
    Vertices(:,1)=Vertices(:,1)/offsets(1);
    Vertices(:,2)=Vertices(:,2)/offsets(2);
    Vertices(:,3)=Vertices(:,3)/offsets(3);
else
    Vertices=Vertices/offsets;
end

Vertices=((invR)*Vertices')';

Vertices(:,1) = Vertices(:,1) - offsetv(1);
Vertices(:,2) = Vertices(:,2) - offsetv(2);
Vertices(:,3) = Vertices(:,3) - offsetv(3);
