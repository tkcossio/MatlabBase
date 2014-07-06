function R=AAM_MakeSearchModel2D(ShapeAppearanceData,ShapeData,AppearanceData,TrainingData,options)


% Structure which will contain all weighted errors of model versus real
% intensities, by several offsets of the parameters
drdp=zeros(size(ShapeAppearanceData.Evectors,2)+4,6,length(TrainingData),length(AppearanceData.g_mean));

% We use the trainingdata images, to train the model. Because we want
% the background information to be included

% Loop through all training images
for i=1:length(TrainingData);
    % Loop through all model parameters, bot the PCA parameters as pose
    % parameters
    for j = 1:size(ShapeAppearanceData.Evectors,2)+4
        if(j<=size(ShapeAppearanceData.Evectors,2))
            % Model parameters, offsets
            de = [-0.5 -0.3 -0.1 0.1 0.3 0.5];
           
            % First we calculate the real ShapeAppearance parameters of the
            % training data set
            c = ShapeAppearanceData.Evectors'*(ShapeAppearanceData.b(:,i) -ShapeAppearanceData.b_mean);
            
            % Standard deviation form the eigenvalue
            c_std = sqrt(ShapeAppearanceData.Evalues(j));
            for k=1:length(de)
                % Offset the ShapeAppearance parameters with a certain
                % value times the std of the eigenvector
                c_offset=c;
                c_offset(j)=c_offset(j)+c_std *de(k);
            
                % Transform back from  ShapeAppearance parameters to Shape parameters  
                b_offset = ShapeAppearanceData.b_mean + ShapeAppearanceData.Evectors*c_offset;
                b1_offset = b_offset(1:(length(ShapeAppearanceData.Ws)));
                b1_offset= inv(ShapeAppearanceData.Ws)*b1_offset;
                x = ShapeData.x_mean + ShapeData.Evectors*b1_offset;
                pos(:,1)=x(1:end/2); 
                pos(:,2)=x(end/2+1:end);


                
                % Transform the Shape back to real image coordinates
                pos=AAM_align_data_inverse2D(pos,TrainingData(i).tform);
                
                % Get the intensities in the real image. Use those
                % intensities to get ShapeAppearance parameters, which
                % are then used to get model intensities
                [g, g_offset]=RealAndModel(TrainingData,i,pos, AppearanceData,ShapeAppearanceData,options,ShapeData);

                % A weighted sum of difference between model an real
                % intensities gives the "intensity / offset" ratio
                w = exp ((-de(k)^2) / (2*c_std^2))/de(k);
                drdp(j,k,i,:)=(g-g_offset)*w;
            end
        else
            % Pose parameters offsets
            j2=j-size(ShapeAppearanceData.Evectors,2);
            switch(j2)
                case 1 % Translation x
                    de = [-2 -1.2 -0.4 0.4 1.2 2]/2;
                case 2 % Translation y
                    de = [-2 -1.2 -0.4 0.4 1.2 2]/2;
                case 3 % Scaling & Rotation Sx
                    de = [-0.2 -.12 -0.04 0.04 0.12 0.2]/2;
                case 4 % Scaling & Rotation Sy
                    de = [-0.2 -.12 -0.04 0.04 0.12 0.2]/2;
            end
            
            for k=1:length(de)
                tform=TrainingData(i).tform;
                switch(j2)
                    case 1 % Translation x
                        tform.offsetv(1)=tform.offsetv(1)+de(k);
                    case 2 % Translation y
                        tform.offsetv(2)=tform.offsetv(2)+de(k);
                    case 3 % Scaling & Rotation Sx
                        tform.offsetsx=tform.offsetsx+de(k);
                    case 4 % Scaling & Rotation Sy
                        tform.offsetsy=tform.offsetsy+de(k);
                end
                
                % From Shape tot real image coordinates, with a certain
                % pose offset
                pos=AAM_align_data_inverse2D(TrainingData(i).CVertices,  tform);
                
                % Get the intensities in the real image. Use those
                % intensities to get ShapeAppearance parameters, which
                % are then used to get model intensities
                [g, g_offset]=RealAndModel(TrainingData,i,pos, AppearanceData,ShapeAppearanceData,options,ShapeData);
             
                % A weighted sum of difference between model an real
                % intensities gives the "intensity / offset" ratio
                w =exp ((-de(k)^2) / (2*2^2))/de(k);
                drdp(j,k,i,:)=(g-g_offset)*w;
            end
        end
    end
end

% Combine the data to the intensity/parameter matrix, 
% using a pseudo inverse
% for i=1:length(TrainingData);
%     drdpt=squeeze(mean(drdp(:,:,i,:),2));
%     R(:,:,i) = (drdpt * drdpt')\drdpt;
% end
% % Combine the data intensity/parameter matrix of all training datasets.
% %
% % In case of only a few images, it will be better to use a weighted mean
% % instead of the normal mean, depending on the probability of the trainingset
% R=mean(R,3);    

drdpt=squeeze(mean(mean(drdp,3),2));
R=pinv(drdpt)';
%R = (drdpt * drdpt')\drdpt;


function [g, g_offset]=RealAndModel(TrainingData,i,pos, AppearanceData,ShapeAppearanceData,options,ShapeData)
   % Sample the image intensities in the training set
   g_offset=AAM_Appearance2Vector2D(TrainingData(i).I,pos, AppearanceData.base_points, AppearanceData.ObjectPixels,ShapeData.TextureSize,ShapeData.Tri);
   g_offset=AAM_NormalizeAppearance2D(g_offset);
   
   % Combine the Shape and Intensity (Appearance) vector 
   b1 = ShapeAppearanceData.Ws * ShapeData.Evectors' * ([TrainingData(i).CVertices(:,1);TrainingData(i).CVertices(:,2)]-ShapeData.x_mean);
   b2 = AppearanceData.Evectors' * (g_offset-AppearanceData.g_mean);
   b = [b1;b2];
   % Calculate the ShapeAppearance parameters
   c2 = ShapeAppearanceData.Evectors'*(b -ShapeAppearanceData.b_mean);
    
   % Go from ShapeAppearance parameters to Appearance parameters
   b = ShapeAppearanceData.b_mean + ShapeAppearanceData.Evectors*c2;
   b2 = b(size(ShapeAppearanceData.Ws,1)+1:end);
   
   % From apperance parameters to intensities
   g = AppearanceData.g_mean + AppearanceData.Evectors*b2;
    
                