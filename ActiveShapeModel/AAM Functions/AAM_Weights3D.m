function Ws = AAM_Weights3D(TrainingData,ShapeData,AppearanceData,options)
nl=length(ShapeData.x_mean)/3;

Change =zeros(length(TrainingData), size(ShapeData.Evectors,2));
for i=1:length(TrainingData)
    % Remove translation and rotation, as done when training the
    % model.
    [pos,tform]=AAM_align_data3D(TrainingData(i).Vertices,ShapeData.MeanVertices,options );
        
    % Describe the model by a vector b with model parameters
    x = [pos(:,1);pos(:,2);pos(:,3)];
    b = ShapeData.Evectors'*(x - ShapeData.x_mean);
    
    % Get the intensities of the untransformed shape.
    % Because noisy eigenvectors from the shape were removed, the 
    % contour is on a little different position and
    % intensities probabbly differ a little bit from the orignal appearance
    x_normal= ShapeData.x_mean + ShapeData.Evectors*b;
    
    pos_normal(:,1) = x_normal(1:nl);
    pos_normal(:,2) = x_normal(nl+1:nl*2);
    pos_normal(:,3) = x_normal(nl*2+1:end);
    
    pos_normal=AAM_align_data_inverse3D(pos_normal,tform,options);
    g_normal=AAM_Appearance2Vector3D(TrainingData(i).I,pos_normal, AppearanceData.base_points, AppearanceData.ObjectPixels,ShapeData.TextureSize,AppearanceData.Tetra,options,ShapeData.Faces);
    g_normal=AAM_NormalizeAppearance3D(g_normal,options);
                     
    for j = 1:size(ShapeData.Evectors,2)
            for k=[-0.5 0.5]
                % Change on model parameter a little bit, to see the influence
                % from the shape parameters on appearance parameters
                b_offset=b;  b_offset(j)=b_offset(j)+k;

                % Transform the model parameter vector b , back to contour positions
                x_offset= ShapeData.x_mean + ShapeData.Evectors*b_offset;
                pos_offset(:,1)=x_offset(1:nl)';
                pos_offset(:,2)=x_offset(nl+1:nl*2)';
                pos_offset(:,3)=x_offset(nl*2+1:end)';
                

                % Now add the previously removed translation and rotation
                [pos_offset]=AAM_align_data_inverse3D(pos_offset,tform,options);

                g_offset=AAM_Appearance2Vector3D(TrainingData(i).I,pos_offset, AppearanceData.base_points, AppearanceData.ObjectPixels,ShapeData.TextureSize,AppearanceData.Tetra,options,ShapeData.Faces);
                g_offset=AAM_NormalizeAppearance3D(g_offset,options);
                Change(i,j) = Change(i,j)+double(sqrt(sum((g_offset-g_normal).^2)/length(g_normal)));
            end
    end
end

Ws=zeros(size(ShapeData.Evectors,2),size(ShapeData.Evectors,2));
for j = 1:size(ShapeData.Evectors,2),
     Ws(j,j) = mean(Change(:,j));
end
