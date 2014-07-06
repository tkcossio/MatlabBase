function Ws = AAM_Weights2D(TrainingData,ShapeData,AppearanceData,options)

Change =zeros(length(TrainingData), size(ShapeData.Evectors,2));
for i=1:length(TrainingData)
    % Remove translation and rotation, as done when training the
    % model.
    [pos,tform]=AAM_align_data2D(TrainingData(i).Vertices,ShapeData.MeanVertices );
        
    % Describe the model by a vector b with model parameters
    x = [pos(:,1);pos(:,2)];
    b = ShapeData.Evectors'*(x - ShapeData.x_mean);
    
    % Get the intensities of the untransformed shape.
    % Because noisy eigenvectors from the shape were removed, the 
    % contour is on a little different position and
    % intensities probabbly differ a little bit from the orignal appearance
    x_normal= ShapeData.x_mean + ShapeData.Evectors*b;
    pos_normal(:,1)=x_normal(1:end/2)'; 
    pos_normal(:,2)=x_normal(end/2+1:end)';
    pos_normal=AAM_align_data_inverse2D(pos_normal,tform);
    g_normal=AAM_Appearance2Vector2D(TrainingData(i).I,pos_normal, AppearanceData.base_points, AppearanceData.ObjectPixels,ShapeData.TextureSize,ShapeData.Tri);
    g_normal=AAM_NormalizeAppearance2D(g_normal);
                     
    for j = 1:size(ShapeData.Evectors,2)
            for k=[-0.5 0.5]
                % Change on model parameter a little bit, to see the influence
                % from the shape parameters on appearance parameters
                b_offset=b;  b_offset(j)=b_offset(j)+k;

                % Transform the model parameter vector b , back to contour positions
                x_offset= ShapeData.x_mean + ShapeData.Evectors*b_offset;
                pos_offset(:,1)=x_offset(1:end/2)';
                pos_offset(:,2)=x_offset(end/2+1:end)';

                % Now add the previously removed translation and rotation
                [pos_offset]=AAM_align_data_inverse2D(pos_offset,tform);

                g_offset=AAM_Appearance2Vector2D(TrainingData(i).I,pos_offset, AppearanceData.base_points, AppearanceData.ObjectPixels,ShapeData.TextureSize,ShapeData.Tri);
                g_offset=AAM_NormalizeAppearance2D(g_offset);
                Change(i,j) = Change(i,j)+sqrt(sum((g_offset-g_normal).^2)/length(g_normal));
            end
    end
end

Ws=zeros(size(ShapeData.Evectors,2),size(ShapeData.Evectors,2));
for j = 1:size(ShapeData.Evectors,2),
     Ws(j,j) = mean(Change(:,j));
end
