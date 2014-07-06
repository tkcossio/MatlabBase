% Load the trained JAW active shape model
load('jaw_aam_data');

% Select and load a Dataset to test the model
i=2;
is=num2str(i); number = '000'; number(end-length(is)+1:end)=is; 
filename=['Images3D\segm' number '.mat'];
load(filename);
Itest=single(V)-0.5*single(imerode(V,ones(5,5,5))); Itest=Itest+rand(size(Itest))*0.1;

% Initial position offset and rotation, of the initial/mean contour
scale=1;
tform= Data{scale}.ShapeData.MeantForm;

% Apply the ASM model onm the test image
[posV,I_model,I_segment]=AAM_ApplyModel3D(Itest,tform,Data,options);

FV.vertices=posV;
FV.faces=Data{scale}.ShapeData.Faces;
showcs3(Itest)
hold on; patch(FV,'facecolor',[0 0 1],'edgecolor', 'none'); axis('vis3d'); view(3); camlight;

showcs3(I_model)

showcs3(Itest/2+I_segment/2);
