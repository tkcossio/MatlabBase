function [greyvector,J]=AAM_Appearance2Vector2D(I,Vertices,base_points,ObjectPixels,texturesize,tri)
% Transform the hands images first into the mean texture image, and than
% transform the image into a vector
%
% greyvector=Appearance2Vector(base_points,Vertices, ObjectPixels,texturesize)
%
%
warning('off','Images:cp2tform:foldOverTriangles');

% The input image coordinates of a training set
input_points = Vertices;

% Make the transformation structure, note that x and y are switched
% because Matlab uses the second dimensions as x and first as y.
% Piecewise-Linear 
xy=[input_points(:,2) input_points(:,1)];
uv=[base_points(:,2) base_points(:,1)];

% Transform the image into the default texture image
if(true)
    % Remove control points which give folded over triangles with cp2tform
    %[xy uv]=PreProcessCp2tform(xy,uv);
    %trans_prj = cp2tform(xy,uv,'piecewise linear');
    %J = imtransform(I,trans_prj,'Xdata',[1 texturesize(1)],'YData',[1 texturesize(2)],'XYscale',1);
    %J(isnan(J))=0;

    J = warp_triangle(I,xy,uv,texturesize,tri);
  
else
    [O_trans,Spacing]=point_registration(texturesize,[uv(:,2) uv(:,1)],[xy(:,2) xy(:,1)],struct('MaxRef',6));
    J = bspline_transform(O_trans,I,Spacing,3,texturesize);
end

  

% Store the transformed texture as a vector
if(size(I,3)==1)
    greyvector=J(ObjectPixels);
else
    Jr=J(:,:,1); Jg=J(:,:,2); Jb=J(:,:,3);
    greyvector=[Jr(ObjectPixels);Jg(ObjectPixels);Jb(ObjectPixels)];
end



