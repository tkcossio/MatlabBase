% This Script shows an example of a working basic Active Appearance Model (AAM),
% with a few hand pictures.
%
% Literature used:
% - T.F. Cootes, G.J Edwards, and C,J. Taylor "Active Appearance Models",
%   Proc. European Conference on Computer Vision 1998
% - T.F. Cootes, G.J Edwards, and C,J. Taylor "Active Appearance Models",
%   IEEE Transactions on Pattern Analysis and Machine Intelligence 2001
%
% Functions are written by D.Kroon University of Twente (March 2010)

%Add functions path to matlab search path
functionname='AAM_2D_example.m'; functiondir=which(functionname);
functiondir=functiondir(1:end-length(functionname));
addpath([functiondir 'AAM Functions'])
addpath([functiondir 'Functions'])
addpath([functiondir 'PieceWiseLinearWarp_version2'])
addpath([functiondir 'PieceWiseLinearWarp_version2/functions'])


% Try to compile c-files
cd([functiondir 'PieceWiseLinearWarp_version2/functions'])
try
    mex('warp_triangle_double.c','image_interpolation.c');
catch ME
    disp('compile c-files failed: example will be slow');
end
cd(functiondir);

%% Set options
% Number of contour points interpolated between the major landmarks.
options.ni=20;
% Set normal appearance/contour, limit to +- m*sqrt( eigenvalue )
options.m=3;
% Size of appearance texture as amount of orignal image
options.texturesize=1;
% If verbose is true all debug images will be shown.
options.verbose=true;
% Number of image scales
options.nscales=4;
% Number of search itterations
options.nsearch=15;

%% Load training data
% First Load the Hand Training DataSets (Contour and Image)
% The LoadDataSetNiceContour, not only reads the contour points, but
% also resamples them to get a nice uniform spacing, between the important
% landmark contour points.
TrainingData=struct;
for i=1:10
    is=num2str(i); number = '000'; number(end-length(is)+1:end)=is;
    filename=['Fotos/train' number '.mat'];
    [TrainingData(i).Vertices, TrainingData(i).Lines]=LoadDataSetNiceContour(filename,options.ni,options.verbose);
    filename=['Fotos/train' number '.jpg'];
    
    I=im2double(imread(filename));
    
    if(options.verbose)
        Vertices=TrainingData(i).Vertices;
        Lines=TrainingData(i).Lines;
        t=mod(i-1,4); if(t==0), figure; end
        subplot(2,2,t+1), imshow(I); hold on;
        P1=Vertices(Lines(:,1),:); P2=Vertices(Lines(:,2),:);
        plot([P1(:,2) P2(:,2)]',[P1(:,1) P2(:,1)]','b');
        drawnow;
    end
    TrainingData(i).I=I;
end

% The Active Appearance Model is constructed for multiple image scales.
% During search first the coarse scale is used to detect the object
% followed by the finer scales (larger images). This makes the model robust
% against initial location and local minima

% The structure which will contain the AAM model for 4 image scales
Data=cell(1,4);
for scale=1:options.nscales
    %% Shape Model %%
    % Make the Shape model, which finds the variations between contours
    % in the training data sets. And makes a PCA model describing normal
    % contours
    [ShapeData,TrainingData] = AAM_MakeShapeModel2D(TrainingData,options);
    
    % Show some eigenvector variations
    if(options.verbose)
        figure,
        for i=1:6
            xtest = ShapeData.x_mean + ShapeData.Evectors(:,i)*sqrt(ShapeData.Evalues(i))*3;
            subplot(2,3,i), hold on;
            plot(xtest(end/2+1:end),xtest(1:end/2),'r.');
            plot(ShapeData.x_mean(end/2+1:end),ShapeData.x_mean(1:end/2),'b.');
        end
    end
    
    %% Appearance model %%
    % Piecewise linear image transformation is used to align all texture
    % information inside the object (hand), to the mean handshape.
    % After transformation of all trainingdata textures to the same shape
    % PCA is used to describe the mean and variances of the object texture.
    AppearanceData=AAM_MakeAppearanceModel2D(TrainingData,ShapeData,options);
    
    % Show some appearance mean and eigenvectors
    if(options.verbose)
        figure,
        I_texture=AAM_Vector2Appearance2D(AppearanceData.g_mean,AppearanceData.ObjectPixels,ShapeData.TextureSize);
        subplot(2,2,1),imshow(I_texture,[]); title('mean grey');
        I_texture=AAM_Vector2Appearance2D(AppearanceData.Evectors(:,1),AppearanceData.ObjectPixels,ShapeData.TextureSize);
        subplot(2,2,2),imshow(I_texture,[]); title('first eigenv');
        I_texture=AAM_Vector2Appearance2D(AppearanceData.Evectors(:,2),AppearanceData.ObjectPixels,ShapeData.TextureSize);
        subplot(2,2,3),imshow(I_texture,[]); title('second eigenv');
        I_texture=AAM_Vector2Appearance2D(AppearanceData.Evectors(:,3),AppearanceData.ObjectPixels,ShapeData.TextureSize);
        subplot(2,2,4),imshow(I_texture,[]); title('third eigenv');
    end
    
    %% Combined model %%
    % Often Shape and Texture are correlated in some way. Thus we can use
    % PCA to get a combined Shape-Appearance model.
    ShapeAppearanceData=AAM_CombineShapeAppearance2D(TrainingData,ShapeData,AppearanceData,options);
    
    if (options.verbose);
        % Show diffference between an original texture, and the original
        % texture described by the combined model. (Are exactly the same
        % if no eigenvectors variations are removed after PCA)
        for i=1:2
            % Orignal texture
            I_texture=AAM_Vector2Appearance2D(AppearanceData.g(:,i),AppearanceData.ObjectPixels,ShapeData.TextureSize);
            
            % Appearance parameters to shape-appearance parameters
            c = ShapeAppearanceData.Evectors'*(ShapeAppearanceData.b(:,i) -ShapeAppearanceData.b_mean) ;
            % Back from Shape-appearance parameters to appearance
            % parameters
            b = ShapeAppearanceData.b_mean + ShapeAppearanceData.Evectors*c;
            b2 = b((length(ShapeAppearanceData.Ws)+1):end);
            g = AppearanceData.g_mean + AppearanceData.Evectors*b2;
            I_texture2=AAM_Vector2Appearance2D(g,AppearanceData.ObjectPixels,ShapeData.TextureSize);
            
            figure,
            subplot(1,3,1),imshow(I_texture,[]); title('original texture');
            subplot(1,3,2),imshow(I_texture2,[]); title('Texture made by the combined model');
            subplot(1,3,3),imshow((I_texture2-I_texture).^2,[]);
            title('Difference between original texture and texture described by the model');
        end
    end
    
    %% Search Model %%
    % The Search Model is used to find the object location and
    % shape-appearance parameters, in a test set.
    % Training is done by displacing the Model and translation parameters
    % with a known amount, and measuring the error, between the intensities
    % form the real image, and those intensities described by the model.
    % The found error correlations are used to make the inverse model which
    % gives the optimial parameter update and new location when you input
    % the error vector with difference between model and real intensities.
    R=AAM_MakeSearchModel2D(ShapeAppearanceData,ShapeData,AppearanceData,TrainingData,options);
    
    % The PCA model is finished, and we store the resulting variables
    % in the data structure, which will contain the Model for 4 different
    % image scales
    Data{scale}.R=R;
    Data{scale}.ShapeAppearanceData=ShapeAppearanceData;
    Data{scale}.ShapeData=ShapeData;
    Data{scale}.AppearanceData=AppearanceData;
    Data{scale}.TrainingData=TrainingData;
    
    % Transform the image to a coarser scale, and update the
    % contour positions.
	if(scale<options.nscales)
		for i=1:length(TrainingData)
			TrainingData(i).Vertices=(TrainingData(i).Vertices-0.5)/2+0.5;
			TrainingData(i).I=imresize(TrainingData(i).I,1/2);
		end
	end
end


%% Applying the model %%%
% Finding and describing the object (hand) in a test image, using
% the multiscale model. We start with the mean shape on the location
% specified below, and use the searchmodel to find the hand location,
% shape and appearance.
Itest=(im2double(imread('Fotos/test001.jpg')));

% Initial position offset and rotation, of the initial/mean contour
tform.offsetsx=1; tform.offsetsy=0; tform.offsetr=0; tform.offsetv=[0 0];
pos(:,1)=Data{1}.ShapeData.x_mean(1:end/2)';
pos(:,2)=Data{1}.ShapeData.x_mean(end/2+1:end)';
% Parameters describing, scale and rotation in comparison with the mean
% hand shape, sx = size*cos(angle), sy=size*sin(angle)

pos=AAM_align_data_inverse2D(pos,tform);

% Select the best starting position with the mouse
[x,y]=SelectPosition(Itest,pos(:,1),pos(:,2));
%x = 190; y = 412;
tform.offsetv=[-x -y];

AAM_ApplyModel2D(Itest,tform,Data,options);



