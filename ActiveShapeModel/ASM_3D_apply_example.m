% Load the trained JAW active shape model
load('jaw_asm_data');

% Select and load a Dataset to test the model
i=2;
is=num2str(i); number = '000'; number(end-length(is)+1:end)=is; 
filename=['Images3D\segm' number '.mat'];
load(filename);
Itest=single(V)-0.5*single(imerode(V,ones(5,5,5))); Itest=Itest+rand(size(Itest))*0.1;

% Initial position offset and rotation, of the initial/mean contour
tform= ShapeData.MeantForm;
% tform.offsetV= [-213.6123 -200.7283 -135.9106]
% tform.offsetrxy= -0.0284
% tform.offsetryz= 1.5486;

% Apply the ASM model onm the test image
posV=ASM_ApplyModel3D(Itest,tform,ShapeData,AppearanceData,options);

% Show the resulting shape
FV.vertices=posV; FV.faces=Faces;
figure, patch(FV,'facecolor',[0 0 1],'edgecolor', 'none'); axis equal; view(3); camlight



