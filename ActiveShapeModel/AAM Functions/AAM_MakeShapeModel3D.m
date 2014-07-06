function [ShapeData TrainingData]= AAM_MakeShapeModel3D(TrainingData,options)
% Number of datasets
s=length(TrainingData);

% Number of landmarks
nl = size(TrainingData(1).Vertices,1);

%% Shape model
MeanVertices=TrainingData(1).Vertices;

% Remove rotation and translation and scale : Procrustes analysis 
AllVertices=zeros([size(MeanVertices) s]);
Alloffsetv=zeros([3 s]);
Alloffsetq=zeros([4 s]);
if(options.scale3)
    Alloffsets=zeros([3 s]);
else
    Alloffsets=zeros([1 s]);
end

for k=1:2
    for i=1:s
        [TrainingData(i).CVertices, TrainingData(i).tform]=AAM_align_data3D(TrainingData(i).Vertices,MeanVertices,options);
        AllVertices(:,:,i)=TrainingData(i).CVertices;
        Alloffsetv(:,i)=TrainingData(i).tform.offsetv;
        Alloffsetq(:,i)=TrainingData(i).tform.offsetq;
        Alloffsets(:,i)=TrainingData(i).tform.offsets;
    end
    tform=struct;
    tform.offsetv=mean(Alloffsetv,2)';
    tform.offsetq=mean(Alloffsetq,2)';
    tform.offsets=mean(Alloffsets,2)';
    CVertices=mean(AllVertices,3);
    MeanVertices=AAM_align_data_inverse3D(CVertices,tform,options);
end
Meantform=tform;
for i=1:s
    [TrainingData(i).CVertices, TrainingData(i).tform]=AAM_align_data3D(TrainingData(i).Vertices,MeanVertices,options);
end

% Construct a matrix with all contour point data of the training data set
x=zeros(nl*3,s);
for i=1:length(TrainingData)
    x(:,i)=[TrainingData(i).CVertices(:,1);TrainingData(i).CVertices(:,2);TrainingData(i).CVertices(:,3)];
end

[Evalues, Evectors, x_mean]=PCA(x);

% Keep only 98% of all eigen vectors, (remove contour noise)
i=find(cumsum(Evalues)>sum(Evalues)*0.99,1,'first'); 
Evectors=Evectors(:,1:i);
Evalues=Evalues(1:i);

% Calculate variances in rotation and scale
q=zeros(s,4);
for i=1:s
    q(i,:)=TrainingData(i).tform.offsetq;
end
varq=var(q,0,1);
varq=options.varq;

% Store the Eigen Vectors and Eigen Values
ShapeData.Evectors=Evectors;
ShapeData.Evalues=Evalues;
ShapeData.x_mean=x_mean;
ShapeData.x = x;
ShapeData.QVariance = varq;
ShapeData.TVariance = [2 2 2];
ShapeData.SVariance = [0.1 0.1  0.1];
ShapeData.MeanVertices = MeanVertices;
ShapeData.Faces = TrainingData(1).Faces;
ShapeData.MeantForm=Meantform;
if(length(options.texturesize)==1)
	ts=ceil(max(max(ShapeData.x_mean(:)),-min(ShapeData.x_mean(:)))*2*options.texturesize);
	ShapeData.TextureSize=[ts ts ts];
else
	ShapeData.TextureSize=options.texturesize;
end



