function AppearanceData=AAM_MakeAppearanceModel3D(TrainingData,ShapeData,options)
% Make the gray-level Appearance Model

% Coordinates of mean contour
nl=length(ShapeData.x_mean)/3;
base_points =[ShapeData.x_mean(1:nl) ShapeData.x_mean(nl+1:nl*2)  ShapeData.x_mean(nl*2+1:end)];

if(isempty(ShapeData.Faces))
	[~,F] = plot3t(base_points(:,1),base_points(:,2),base_points(:,3),6,'r',8); 
	base_points=F.vertices;
	Faces=F.faces;
else
    Faces=ShapeData.Faces;
end
xyz=base_points;


% Normalize the base points to range 0..1
base_points = base_points - repmat(min(base_points),size(base_points,1),1);
base_points = base_points ./ repmat(max(base_points),size(base_points,1),1);

% Transform the mean contour points into the coordinates in the texture
% image.
base_points(:,1)=1+options.borderpixel+(ShapeData.TextureSize(1)-1-2*options.borderpixel)*base_points(:,1);
base_points(:,2)=1+options.borderpixel+(ShapeData.TextureSize(2)-1-2*options.borderpixel)*base_points(:,2);
base_points(:,3)=1+options.borderpixel+(ShapeData.TextureSize(3)-1-2*options.borderpixel)*base_points(:,3);
uvw=base_points;

% Draw the surface as one closed white surface and fill the resulting
% (jaw) object

Nuvw=-patchnormals(struct('vertices',base_points,'faces',Faces))*options.borderpixel;

VB=base_points+Nuvw;
ObjectPixels = polygon2voxel(struct('vertices',VB,'faces',Faces),ShapeData.TextureSize,'none', false);    
ObjectPixels = imfill(ObjectPixels ,'holes');

% Calculate triangulation (tetrahedrons), including border 
if(options.constanttetra)
	FVxyz.vertices=xyz; FVxyz.faces=Faces;
    scale=mean(sqrt(sum(bsxfun(@minus,xyz,mean(xyz,2)).^2,2)))/mean(sqrt(sum(bsxfun(@minus,uvw,mean(uvw,2)).^2,2)));
    Nxyz=-patchnormals(FVxyz)*(options.borderpixel+1)*scale;
	Tetra= delaunayn([xyz;xyz+Nxyz],{'Qt','Qbb','Qc','Qz'});
else
	Tetra=[];
end

% Number of datasets
s=length(TrainingData);

% Transform the hands images first into the mean texture image, and than
% transform the image into a vector using the pixellocations of the object
% found here above.

% Construct a matrix with all appearance data of the training data set
npixels=sum(ObjectPixels(:)>0);
g=zeros(npixels,s);
for i=1:s
    g(:,i)=AAM_Appearance2Vector3D(TrainingData(i).I,TrainingData(i).Vertices, base_points, ObjectPixels,ShapeData.TextureSize,Tetra,options,ShapeData.Faces);
end

% Normalize the greylevels, to compensate for illumination 
for i=1:s
    g(:,i)=AAM_NormalizeAppearance3D(g(:,i),options);
end
[Evalues, Evectors, g_mean]=PCA(g);

% Keep only 99% of all eigen vectors, (remove contour noise)
i=find(cumsum(Evalues)>sum(Evalues)*0.99,1,'first'); 
Evectors=Evectors(:,1:i);
Evalues=Evalues(1:i);

% Store the Eigen Vectors and Eigen Values
AppearanceData.Evectors=Evectors;
AppearanceData.Evalues=Evalues;
AppearanceData.g_mean=g_mean;
AppearanceData.g = g;
AppearanceData.ObjectPixels=ObjectPixels;
AppearanceData.base_points=base_points;
AppearanceData.Tetra=Tetra;
