% Train an active Shape Model of the JAW
%
% For speed-up compile the function patchnormals_double.c in ASM functions
% using  mex patchnormals_double.c -v
%
% Functions are written by D.Kroon University of Twente (March 2011)

% Add functions path to matlab search path
functionname='ASM_3D_train_example.m'; functiondir=which(functionname);
functiondir=functiondir(1:end-length(functionname));
addpath([functiondir 'Functions'])
addpath([functiondir 'ASM Functions'])
addpath([functiondir 'polygon2voxel_version1j'])
addpath([functiondir 'PatchNormals_version1'])

% Compile c-files
cd([functiondir 'InterpFast_version1'])
mex('interp3fast_double.c','image_interpolation.c');
mex('interp3fast_single.c','image_interpolation.c');
cd([functiondir 'polygon2voxel_version1j'])
mex('polygon2voxel_double.c');
cd([functiondir 'PatchNormals_version1'])
mex('patchnormals_double.c');
cd(functiondir);





%% Set options
% If verbose is true all debug images will be shown.
% Length of landmark intensity profile
options.k = 8; 
% Search length (in pixels) for optimal contourpoint position, 
% in both normal directions of the contourpoint.
options.ns=6;
% Number of image resolution scales
options.nscales=3;
% Set normal contour, limit to +- m*sqrt( eigenvalue )
options.m=3;
% Number of search itterations
options.nsearch=[5 5 5];
% If verbose is true all debug images will be shown.
options.verbose=true;
% The original minimal Mahanobis distance using edge gradient (true)
% or new minimal PCA parameters us
options.originalsearch=false;
% During search try multiple initial positions
options.optimizestart=false;

%% Load training data
if(options.verbose), disp('loading training data'); drawnow; end

TrainingData=struct;
for i=1:10
    is=num2str(i); number = '000'; number(end-length(is)+1:end)=is; 
    filename=['Images3D\surface' number '.mat'];
    load(filename);
    filename=['Images3D\segm' number '.mat'];
    load(filename);
    % Make a fake dataset of the JAW
    I=single(V)-0.5*single(imerode(V,ones(5,5,5))); I=I+rand(size(I))*0.1;
    TrainingData(i).Vertices=Vertices;
    TrainingData(i).I=I;
end
TrainingData(1).Faces=Faces;
clear I; clear V; clear Vertices; clear RefSurface;

if(options.verbose),
    FV.vertices=TrainingData(2).Vertices; FV.faces=Faces;
    showcs3(TrainingData(2).I); hold on;
    patch(FV,'facecolor',[1 0 0],'facealpha',0.5,'edgecolor','none')
end

%% Shape Model %%
% Make the Shape model, which finds the variations between contours
% in the training data sets. And makes a PCA model describing normal
% contours
if(options.verbose), disp('Creating Shape Model'); drawnow; end
[ShapeData TrainingData]= ASM_MakeShapeModel3D(TrainingData);

if(options.verbose),
    nl=length(ShapeData.x_mean)/3;
    posV(:,1)=ShapeData.x_mean(1:nl)'; 
    posV(:,2)=ShapeData.x_mean(nl+1:nl*2)';
    posV(:,3)=ShapeData.x_mean(nl*2+1:end)';
    FV.vertices=posV;
    FV.faces=Faces;
    figure, patch(FV,'facecolor',[0 0 1],'edgecolor', 'none'); axis('vis3d'); view(3); camlight
end


%% Appearance model %%
% Make the Appearance model, which samples a intensity pixel profile/line 
% perpendicular to each contourpoint in each trainingdataset. Which is 
% used to build correlation matrices for each landmark. Which are used
% in the optimization step, to find the best fit.
if(options.verbose), disp('Creating Appearance model'); drawnow; end
AppearanceData = ASM_MakeAppearanceModel3D(TrainingData,Faces,options);

if(options.verbose), disp('Storing Data'); drawnow; end
save('jaw_asm_data','-v7.3','ShapeData','AppearanceData','options','Faces');
