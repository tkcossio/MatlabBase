function ShapeAppearanceData=AAM_CombineShapeAppearance3D(TrainingData,ShapeData,AppearanceData,options)
% This functions combines the shape and appearance of the objects, by
% adding the weighted vector describing shape and appearance, followed by
% PCA

% Get weight matrix. The Weights are a scaling between texture and shape
% to give a change in shape parameters approximately the same 
% influences as texture parameters.
Ws = AAM_Weights3D(TrainingData,ShapeData,AppearanceData,options);

% Combine the Contour and Appearance data
b=zeros(size(ShapeData.Evectors,2)+size(AppearanceData.Evectors,2),length(TrainingData));
for i=1:length(TrainingData)
    b1 = Ws * ShapeData.Evectors' * (ShapeData.x(:,i)-ShapeData.x_mean);
    b2 = AppearanceData.Evectors' * (AppearanceData.g(:,i)-AppearanceData.g_mean);
    b(:,i)=[b1;b2];
end


% By definition b_mean is zero, because the shapedata en appearancedata
% already had the mean substracted
[Evalues, Evectors, b_mean]=PCA(b);

% Keep only 99% of all eigen vectors, (remove noise)
i=find(cumsum(Evalues)>sum(Evalues)*0.99,1,'first'); 
Evectors=Evectors(:,1:i);
Evalues=Evalues(1:i);

ShapeAppearanceData.Evectors=Evectors;
ShapeAppearanceData.Evalues=Evalues;
ShapeAppearanceData.b_mean=b_mean;
ShapeAppearanceData.b=b;
ShapeAppearanceData.Ws=Ws;
