function [ShapeData TrainingData]= ASM_MakeShapeModel3D(TrainingData)

% Number of datasets
s=length(TrainingData);

% Number of landmarks
nl = size(TrainingData(1).Vertices,1);

%% Shape model

% Remove rotation and translation 
% (Procrustes analysis would also be possible, see AAM_align_data)
MeanVertices=TrainingData(1).Vertices;
offsetryz=0; offsetrxy=0; offsetV=[0 0 0];
for i=1:s
    [TrainingData(i).CVertices, TrainingData(i).tform]=ASM_align_data3D(TrainingData(i).Vertices,MeanVertices);
    offsetV=offsetV+TrainingData(i).tform.offsetV;
    offsetrxy=offsetrxy+TrainingData(i).tform.offsetrxy;
    offsetryz=offsetryz+TrainingData(i).tform.offsetryz;
end
MeantForm.offsetV=offsetV/s;
MeantForm.offsetrxy=offsetrxy/s;
MeantForm.offsetryz=offsetryz/s;

% Construct a matrix with all contour point data of the training data set
x=zeros(nl*3,s);
for i=1:length(TrainingData)
    x(:,i)=[TrainingData(i).CVertices(:,1);TrainingData(i).CVertices(:,2);TrainingData(i).CVertices(:,3)];
end

[Evalues, Evectors, x_mean]=PCA(x);

% Keep only 98% of all eigen vectors, (remove contour noise)
i=find(cumsum(Evalues)>sum(Evalues)*0.98,1,'first'); 
Evectors=Evectors(:,1:i);
Evalues=Evalues(1:i);

% Store the Eigen Vectors and Eigen Values
ShapeData.Evectors=Evectors;
ShapeData.Evalues=Evalues;
ShapeData.x_mean=x_mean;
ShapeData.x=x;
ShapeData.MeanVertices = MeanVertices;
ShapeData.MeantForm =MeantForm;
ShapeData.Faces=TrainingData(1).Faces;




