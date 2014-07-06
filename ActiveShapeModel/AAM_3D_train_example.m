% Train an Active Appearance Model of the JAW
%
% Functions are written by D.Kroon University of Twente (March 2011)

% Add functions path to matlab search path
functionname='AAM_3D_train_example.m'; functiondir=which(functionname);
functiondir=functiondir(1:end-length(functionname));
addpath([functiondir 'AAM Functions'])
addpath([functiondir 'Functions'])
addpath([functiondir 'PieceWiseLinearWarp_version2'])
addpath([functiondir 'PieceWiseLinearWarp_version2\functions'])
addpath([functiondir 'polygon2voxel_version1j'])
addpath([functiondir 'PatchNormals_version1'])
addpath([functiondir 'InterpFast_version1'])


% Compile c-files
cd([functiondir 'PieceWiseLinearWarp_version2'])
compile_c_files
cd([functiondir 'InterpFast_version1'])
compile_c_files

cd([functiondir 'polygon2voxel_version1j'])
mex('polygon2voxel_double.c');
cd(functiondir);
cd([functiondir 'PatchNormals_version1'])
mex('patchnormals_double.c');
cd(functiondir);

%% Set options
% Number of contour points interpolated between the major landmarks.
options.ni=20;
% Set normal appearance/contour, limit to +- m*sqrt( eigenvalue )
options.m=3;
% Size of texture appereance image
%options.texturesize=1;
options.texturesize=[150 150 150];
% If verbose is true all debug images will be shown.
options.verbose=true;
% Number of image scales
options.nscales=4;
% Number of search itterations
options.nsearch=15;
% Scale Q-matrix
%options.scaleqmatrix=10;
options.scaleqmatrix=10;
% Use one time delaunay for all transformations.
% options.constanttetra=true;
options.constanttetra=true;
% Use pose variance in search model building
%options.posevariance=true;
options.posevariance=true;
% Use pinv or orignal pseudo-inverse
%options.usepinv=true;
options.usepinv=true;
% Calculate R for every dataset or collect first mean
%options.allr=false;
options.allr=false;
% Use model parameters of search or directly from image
options.usemodelc=false;
% Optimize c-parameters before starting or after ending
options.optimizecstart=true;
options.optimizecend=true;
options.optimizecmiddle=false;
options.borderpixel=5;
% Normalize grey values
options.normalizeg=true;
% Variation in quaternion vector
options.varq=[0.2 0.2 0.2 0.2];
options.scale3=false;
options.weightShapeAppearance=1;
options.boostnervesigma=10;
options.boostnervesoffset=0.1;
% Use model parameters of search or directly from image
options.usemodelc=true;
options.logerror=false;
options.sigmas=0.1;
options.optimizestart=false;

%% Load training data
if(options.verbose), disp('loading training data'); drawnow; end

TrainingData=struct;
for i=1:5
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
clear I; clear V; clear Vertices; clear RefSurface;
TrainingData(1).Faces=Faces;

if(options.verbose),
    FV.vertices=TrainingData(2).Vertices; FV.faces=Faces;
    showcs3(TrainingData(2).I); hold on;
    patch(FV,'facecolor',[1 0 0],'facealpha',0.5,'edgecolor','none')
end

%% Shape Model %%
% Make the Shape model, which finds the variations between contours
% in the training data sets. And makes a PCA model describing normal
% contours

% The structure which will contain the AAM model for 4 image scales
Data=cell(1,4);
for scale=1:options.nscales
	if(options.verbose), disp(['Scale ' num2str(scale) ' of ' num2str(options.nscales)]); drawnow; end
	
     %% Shape Model %%
    % Make the Shape model, which finds the variations between contours
    % in the training data sets. And makes a PCA model describing normal
    % contours
    if(options.verbose), disp('Creating Shape Model'); drawnow; end
    [ShapeData,TrainingData] = AAM_MakeShapeModel3D(TrainingData,options);
    
    if(options.verbose),
        nl=length(ShapeData.x_mean)/3;
        posV(:,1)=ShapeData.x_mean(1:nl)'; 
        posV(:,2)=ShapeData.x_mean(nl+1:nl*2)';
        posV(:,3)=ShapeData.x_mean(nl*2+1:end)';
        posV=AAM_align_data_inverse3D(posV,TrainingData(2).tform,options);
        
        FV.vertices=posV;
        FV.faces=Faces;
        showcs3(TrainingData(2).I)
        hold on; patch(FV,'facecolor',[0 0 1],'edgecolor', 'none'); axis('vis3d'); view(3); camlight;
        pause(1);
    end
        
    %% Appearance model %%
    % Piecewise linear image transformation is used to align all texture
    % information inside the object (hand), to the mean handshape.
    % After transformation of all trainingdata textures to the same shape
    % PCA is used to describe the mean and variances of the object texture.
	if(options.verbose), disp('Creating Appearance Model'); drawnow; end
    AppearanceData=AAM_MakeAppearanceModel3D(TrainingData,ShapeData,options);
    
  %% Combined model %%
    % Often Shape and Texture are correlated in some way. Thus we can use
    % PCA to get a combined Shape-Appearance model.
	if(options.verbose), disp('Combine Shape and Appearance Model'); drawnow; end
    ShapeAppearanceData=AAM_CombineShapeAppearance3D(TrainingData,ShapeData,AppearanceData,options);
 
  %% Search Model %%
    % The Search Model is used to find the object location and
    % shape-appearance parameters, in a test set.
    % Training is done by displacing the Model and translation parameters
    % with a known amount, and measuring the error, between the intensities
    % form the real image, and those intensities described by the model.
    % The found error correlations are used to make the inverse model which
    % gives the optimial parameter update and new location when you input
    % the error vector with difference between model and real intensities.
	if(options.verbose), disp('Make search Model'); drawnow; end
    R=AAM_MakeSearchModel3D(ShapeAppearanceData,ShapeData,AppearanceData,TrainingData,options);
    
    % The PCA model is finished, and we store the resulting variables
    % in the data structure, which will contain the Model for 4 different
    % image scales
    Data{scale}.R=R;
    Data{scale}.ShapeAppearanceData=ShapeAppearanceData;
    Data{scale}.ShapeData=ShapeData;
    Data{scale}.AppearanceData=AppearanceData;
    
    % Transform the image to a coarser scale, and update the 
    % contour positions.
	if(scale<options.nscales)
		if(options.verbose), disp('Resize training volumes'); drawnow; end
		for i=1:length(TrainingData)
			TrainingData(i).Vertices=(TrainingData(i).Vertices-0.5)/2+0.5;
			TrainingData(i).I=imresize3d(TrainingData(i).I,1/2,[],'cubic','replicate');
		end
	end
end

if(options.verbose), disp('Storing Data'); drawnow; end
save('jaw_aam_data','-v7.3','Data','options');
