function [uv xy]=PreProcessCp2tform(uv,xy)
% When CP2TFORM is used in piecewise linear mode, it uses triangulation
% to calculate the local image transformation. Sometimes a good 
% triangulation of basepoints gives folded triangles when used
% with the input points. The function CP2TFORM try's to remove 
% the controlpoints causing the folding, but often fails with the error:
%    Eliminated * control point pair(s).
%    Fold-over triangles remain. See CP2TFORM reference page.
%
% This function PreProcessCp2tform removes all control points which cause
% fold-over triangles with cp2tform. 
%
% [input_points,base_points]=PreProcessCp2tform(input_points,base_points);
%
% Function is written by D.Kroon University of Twente (March 2010)

% See if there are triangles after triangulation, which are upside down
tri = delaunay(xy(:,1),xy(:,2));
fold_triangles = FindUpsideDown(xy,uv,tri);

while(~isempty(fold_triangles))
    % Find the vertices in the fold-over triangle list which have the
    % largest angle, most probably causing the folding.
    tri_fold=tri(fold_triangles,:); tri_foldt=tri_fold';
    
    % Convert the folded triangle list, to list of vertex coordinates
    x = xy(:,1); y = xy(:,2);
    vx = reshape(x(tri_fold),size(tri_fold,1),3);
    vy = reshape(y(tri_fold),size(tri_fold,1),3);
    
    % Find the vertices with the largest angle in there triangles
    angle=[abs(dot([vx(:,1)-vx(:,3) vy(:,1)-vy(:,3)]',[vx(:,1)-vx(:,2) vy(:,1)-vy(:,2)]')) 
           abs(dot([vx(:,2)-vx(:,3) vy(:,2)-vy(:,3)]',[vx(:,2)-vx(:,1) vy(:,2)-vy(:,1)]')) 
           abs(dot([vx(:,3)-vx(:,1) vy(:,3)-vy(:,1)]',[vx(:,3)-vx(:,2) vy(:,3)-vy(:,2)]'))];
    [dum index] = min(angle); 
    bad_vertices=tri_foldt(index+(0:length(fold_triangles)-1)*3);
  
    % eliminate vertices causing folding
    bad_vertices = unique(bad_vertices);
    xy(bad_vertices,:) = []; uv(bad_vertices,:) = [];

    % See if there are triangles after triangulation, which are upside down
    tri = delaunay(xy(:,1),xy(:,2));
    fold_triangles = FindUpsideDown(xy,uv,tri);
end

function index = FindUpsideDown(xy,uv,tri)
% look for triangles which are upside down 
x = xy(:,1); y = xy(:,2);
u = uv(:,1); v = uv(:,2);

xx = reshape(x(tri),size(tri,1),3);
yy = reshape(y(tri),size(tri,1),3);
uu = reshape(u(tri),size(tri,1),3);
vv = reshape(v(tri),size(tri,1),3);

% Calculate the normal of the polygon
normalxy=cross([xx(:,2)-xx(:,1) yy(:,2)-yy(:,1) zeros(size(xx,1),1)],[xx(:,3)-xx(:,1) yy(:,3)-yy(:,1) zeros(size(xx,1),1)]);
normaluv=cross([uu(:,2)-uu(:,1) vv(:,2)-vv(:,1) zeros(size(xx,1),1)],[uu(:,3)-uu(:,1) vv(:,3)-vv(:,1) zeros(size(xx,1),1)]);

% Find the vertices with normals pointing in different directions
index = find(normalxy(:,3).*normaluv(:,3)<0);
