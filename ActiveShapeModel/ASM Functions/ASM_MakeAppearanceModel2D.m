%% Gray-Level Appearance Model
function AppearanceData=ASM_MakeAppearanceModel2D(TrainingData,options)

% Number of TrainingData sets
s=length(TrainingData);

% Number of landmarks
nl = size(TrainingData(1).Vertices,1);

% Calculate the normals of the contours of all training data
for i=1:s
  [TrainingData(i).Normals]=ASM_GetContourNormals2D(TrainingData(i).Vertices,TrainingData(i).Lines);
end

% Inverse of covariance matrix sometimes badly scaled
warning('off','MATLAB:nearlySingularMatrix');

AppearanceData=struct;
% Get the landmark profiles for 3 image scales (for multiscale ASM)
for itt_res=1:options.nscales
    scale=1 / (2^(itt_res-1));

    % Get the pixel profiles of every landmark perpendicular to the contour
    for i=1:s
        P = (TrainingData(i).Vertices-0.5)*scale+0.5;
        %P = TrainingData(i).Vertices*scale;
        Ismall=imresize(TrainingData(i).I,scale);
        N = TrainingData(i).Normals;
        [TrainingData(i).GrayProfiles,TrainingData(i).DerivativeGrayProfiles]=ASM_getProfileAndDerivatives2D(Ismall,P,N,options.k);
    end

%     if(options.verbose),
%        disp('Process Intensity profiles');
%        figure, hold on;
%        for i=1:s
%            plot(mean(TrainingData(i).GrayProfiles,2),'Color',rand(1,3));
%        end
%        legend('1','2','3','4','5','6')
%        title('Mean intensity profiles');
%        drawnow;
%     end
    
    % Profile length * space for rgb
    pl=(options.k*2+1)*size(TrainingData(1).I,3);
    
    % Calculate a covariance matrix for all landmarks
    PCAData=struct;
    for j=1:nl
        %% The orginal search method using Mahanobis distance with
        % edge gradient information
        dg=zeros(pl,s);
        for i=1:s, dg(:,i)=TrainingData(i).DerivativeGrayProfiles(:,j); end
        dg_mean=mean(dg,2);
        dg=dg-repmat(dg_mean,1,s);
        % Calculate the coveriance matrix and its inverse
        AppearanceData(itt_res).Landmarks(j).S = cov(dg');
        AppearanceData(itt_res).Landmarks(j).Sinv = inv(AppearanceData(itt_res).Landmarks(j).S);
        AppearanceData(itt_res).Landmarks(j).dg_mean = dg_mean;
        
        %% The new search method using PCA on intensities, and minimizing
        % parameters / the distance to the mean during the search.
        % Make a matrix with all intensities
        g=zeros(pl,s);
        for i=1:s, g(:,i)=TrainingData(i).GrayProfiles(:,j); end
                
        [Evalues, Evectors, Emean]=PCA(g);
        % Keep only 98% of all eigen vectors, (remove contour noise)
		i=find(cumsum(Evalues)>sum(Evalues)*0.98,1,'first'); 
        PCAData(j).Evectors= Evectors(:,1:i);
        PCAData(j).Evalues = Evalues(1:i);
        PCAData(j).Emean = Emean;
        
    end
    AppearanceData(itt_res).PCAData=PCAData;
end


