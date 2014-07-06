function [ShapeData TrainingData]= ASM_MakeShapeModel2D(TrainingData)

% Number of datasets
s=length(TrainingData);

% Number of landmarks
nl = size(TrainingData(1).Vertices,1);

%% Shape model

% Remove rotation and translation 
% (Procrustes analysis would also be possible, see AAM_align_data)
for i=1:s
    [TrainingData(i).CVertices, TrainingData(i).tform]=ASM_align_data2D(TrainingData(i).Vertices);
end

% Construct a matrix with all contour point data of the training data set
x=zeros(nl*2,s);
for i=1:length(TrainingData)
    x(:,i)=[TrainingData(i).CVertices(:,1)' TrainingData(i).CVertices(:,2)']';
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
ShapeData.Lines=TrainingData(1).Lines;



