function greyvector=AAM_Appearance2Vector3D(I,Vertices,base_points,ObjectPixels,texturesize,tetra,options,Faces)
% Transform the hands images first into the mean texture image, and than
% transform the image into a vector
%
% greyvector=Appearance2Vector(base_points,Vertices, ObjectPixels,texturesize)
%
%

    
% The input image coordinates of a training set
input_points = Vertices;

% Make the transformation structure, note that x and y are switched
% because Matlab uses the second dimensions as x and first as y.
% Piecewise-Linear 
if(isempty(Faces))
	[~,F] = plot3t(input_points(:,1),input_points(:,2),input_points(:,3),6,'r',8); 
	input_points=F.vertices;
	Faces=F.faces;
end

xyz=[input_points(:,2) input_points(:,1) input_points(:,3)];
uvw=[base_points(:,2) base_points(:,1) base_points(:,3)];

% Transform the image into the default texture image
if(true)
    FVxyz.vertices=xyz; FVxyz.faces=Faces;
    FVuvw.vertices=uvw; FVuvw.faces=Faces;
    
    scale=mean(sqrt(sum(bsxfun(@minus,xyz,mean(xyz,2)).^2,2)))/mean(sqrt(sum(bsxfun(@minus,uvw,mean(uvw,2)).^2,2)));
    
    Nxyz=patchnormals(FVxyz)*(options.borderpixel+1)*scale;
    Nuvw=patchnormals(FVuvw)*(options.borderpixel+1);
    
	if(options.constanttetra)
		J = warp_tetrahedron(I,[xyz;xyz+Nxyz],[uvw;uvw+Nuvw],texturesize,tetra);
	else
		J = warp_tetrahedron(I,[xyz;xyz+Nxyz],[uvw;uvw+Nuvw],texturesize);
	end
else
    [O_trans,Spacing]=point_registration(texturesize,[uvw(:,2) uvw(:,1) uvw(:,3)],[xyz(:,2) xyz(:,1) xyz(:,3)],struct('MaxRef',6));
    J = bspline_transform(O_trans,I,Spacing,3,texturesize);
end

% Store the transformed texture as a vector
check=ObjectPixels>0;
greyvector=J(check).*ObjectPixels(check);
