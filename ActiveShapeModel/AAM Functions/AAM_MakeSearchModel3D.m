function R=AAM_MakeSearchModel3D(ShapeAppearanceData,ShapeData,AppearanceData,TrainingData,options)
nl=length(ShapeData.x_mean)/3;

if(options.scale3), posen=10; else posen=7; end

% Structure which will contain all weighted errors of model versus real
% intensities, by several offsets of the parameters
drdp=zeros(size(ShapeAppearanceData.Evectors,2)+posen,6,length(AppearanceData.g_mean));

% We use the trainingdata images, to train the model. Because we want
% the background information to be included

		
% Loop through all training images
for i=1:length(TrainingData);
    if(options.verbose),
        disp(['Warping Training Image : ' num2str(i)]);
    end
    % Loop through all model parameters, bot the PCA parameters as pose
    % parameters
    for j = 1:size(ShapeAppearanceData.Evectors,2)+posen;
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
                pos(:,1)=x(1:nl);
                pos(:,2)=x(nl+1:nl*2);
                pos(:,3)=x(nl*2+1:end);
                
                % Transform the Shape back to real image coordinates
                pos=AAM_align_data_inverse3D(pos,TrainingData(i).tform,options);
                
                % Get the intensities in the real image. Use those
                % intensities to get ShapeAppearance parameters, which
                % are then used to get model intensities
                [g, g_offset]=RealAndModel(TrainingData,i,pos, AppearanceData,ShapeAppearanceData,options,ShapeData);
                
                % A weighted sum of difference between model an real
                % intensities gives the "intensity / offset" ratio
                w = exp ((-de(k)^2) / (2*c_std^2))/de(k);
                drdp(j,k,:)=(g-g_offset)*w;
            end
        else
			if(options.scale3)
				p_std = [ShapeData.TVariance(:);ShapeData.QVariance(:);ShapeData.SVariance(:)];
            else
				p_std = [ShapeData.TVariance(:);ShapeData.QVariance(:)];
			end
			
            % Pose parameters offsets
            j2=j-size(ShapeAppearanceData.Evectors,2);
            if(~options.posevariance)
                p_std=[2 2 2 1/30 1/30 1/30 1/30];
                switch(j2)
                    case 1 % Translation x
                        de = [-2 -1.2 -0.4 0.4 1.2 2]/2;
                    case 2 % Translation y
                        de = [-2 -1.2 -0.4 0.4 1.2 2]/2;
                    case 3 % Translation z
                        de = [-2 -1.2 -0.4 0.4 1.2 2]/2;
                    case 4 % Scaling & Rotation
                        de = [-0.2 -.12 -0.04 0.04 0.12 0.2]/2;
                    case 5 % Scaling & Rotation
                        de = [-0.2 -.12 -0.04 0.04 0.12 0.2]/2;
                    case 6 % Scaling & Rotation
                        de = [-0.2 -.12 -0.04 0.04 0.12 0.2]/2;
                    case 7 % Scaling & Rotation
                        de = [-0.2 -.12 -0.04 0.04 0.12 0.2]/2;
					case 8 % Scaling & Rotation
                        de = [-0.2 -.12 -0.04 0.04 0.12 0.2]/2;
                    case 9 % Scaling & Rotation
                        de = [-0.2 -.12 -0.04 0.04 0.12 0.2]/2;
                    case 10 % Scaling & Rotation
                        de = [-0.2 -.12 -0.04 0.04 0.12 0.2]/2;
                end
            end
            for k=1:length(de)
                tform=TrainingData(i).tform;
                switch(j2)
                    case 1 % Translation x
                        tform.offsetv(1)=tform.offsetv(1)+de(k)*p_std(j2);
                    case 2 % Translation y
                        tform.offsetv(2)=tform.offsetv(2)+de(k)*p_std(j2);
                    case 3 % Translation z
                        tform.offsetv(3)=tform.offsetv(3)+de(k)*p_std(j2);
                    case 4 % (Scaling) & Rotation (Quaternion)
                        tform.offsetq(1)=tform.offsetq(1)+de(k)*p_std(j2);
                    case 5 % (Scaling) & Rotation (Quaternion)
                        tform.offsetq(2)=tform.offsetq(2)+de(k)*p_std(j2);
                    case 6 % (Scaling) & Rotation (Quaternion)
                        tform.offsetq(3)=tform.offsetq(3)+de(k)*p_std(j2);
                    case 7 % (Scaling) & Rotation (Quaternion)
                        tform.offsetq(4)=tform.offsetq(4)+de(k)*p_std(j2);
				    case 8  %  Scaling x
                        tform.offsets(1)=tform.offsets(1)+de(k)*p_std(j2);
                    case 9  % Scaling y 
                        tform.offsets(2)=tform.offsets(2)+de(k)*p_std(j2);
                    case 10 % Scaling z
                        tform.offsets(3)=tform.offsets(3)+de(k)*p_std(j2);
                end
                
                % From Shape tot real image coordinates, with a certain
                % pose offset
                pos=AAM_align_data_inverse3D(TrainingData(i).CVertices,  tform, options);
                
                % Get the intensities in the real image. Use those
                % intensities to get ShapeAppearance parameters, which
                % are then used to get model intensities
                [g, g_offset]=RealAndModel(TrainingData,i,pos, AppearanceData,ShapeAppearanceData,options,ShapeData);
                
                % A weighted sum of difference between model an real
                % intensities gives the "intensity / offset" ratio
                if(options.posevariance)
                    w = exp ((-de(k)^2) / (2*p_std(j2)^2))/de(k);
                else
                    w =exp ((-de(k)^2) / (2*2^2))/de(k);
                end
                
                drdp(j,k,:)=(g-g_offset)*w;
            end
        end
    end
    
    % Combine the data to the intensity/parameter matrix,
    % using a pseudo inverse
    
    %     dis=zeros([size(drdpt,1) size(drdpt,1)]);
    %     for ix=1:size(drdpt,1)
    %         for iy=1:size(drdpt,1)
    %             dis(ix,iy)=sum(sqrt((drdpt(ix,:)-drdpt(iy,:)).^2));
    %         end
    %     end
    %     dis
    %     figure, imshow(dis)
    if(options.allr)
        drdpt=squeeze(mean(drdp,2));
        if(options.usepinv)
            R1=pinv(drdpt)';
        else
            R1 = (drdpt * drdpt')\drdpt;
        end
        if(i==1), R = R1; else R = R + R1; end
    else
        if(i==1),
            drdpt=drdp;
        else
            drdpt=drdpt+drdp;
        end
    end
end
if(options.allr)
    % Normalize the inverse intensity/parameter matrix,
    R=R/length(TrainingData);
else
    drdpt=squeeze(mean(drdpt,2));
    drdpt=drdpt/length(TrainingData);
    if(options.usepinv)
        R=pinv(drdpt)';
    else
        R = (drdpt * drdpt')\drdpt;
    end
end

if(sum(isnan(R(:)))>0), keyboard; end

% Combine the data intensity/parameter matrix of all training datasets.
%
% In case of only a few images, it will be better to use a weighted mean
% instead of the normal mean, depending on the probability of the trainingset


function [g, g_offset]=RealAndModel(TrainingData,i,pos, AppearanceData,ShapeAppearanceData,options,ShapeData)
% Sample the image intensities in the training set
g_offset=AAM_Appearance2Vector3D(TrainingData(i).I,pos, AppearanceData.base_points, AppearanceData.ObjectPixels,ShapeData.TextureSize,AppearanceData.Tetra,options,ShapeData.Faces);
g_offset=AAM_NormalizeAppearance3D(g_offset,options);

% Combine the Shape and Intensity (Appearance) vector
b1 = ShapeAppearanceData.Ws * ShapeData.Evectors' * ([TrainingData(i).CVertices(:,1);TrainingData(i).CVertices(:,2);TrainingData(i).CVertices(:,3)]-ShapeData.x_mean);
b2 = AppearanceData.Evectors' * (g_offset-AppearanceData.g_mean); b2=double(b2);
b = [b1;b2];
% Calculate the ShapeAppearance parameters
c2 = ShapeAppearanceData.Evectors'*(b -ShapeAppearanceData.b_mean);

% Go from ShapeAppearance parameters to Appearance parameters
b = ShapeAppearanceData.b_mean + ShapeAppearanceData.Evectors*c2;
b2 = b(size(ShapeAppearanceData.Ws,1)+1:end);

% From apperance parameters to intensities
g = AppearanceData.g_mean + AppearanceData.Evectors*b2;

