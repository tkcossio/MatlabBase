function J=warp_triangle_double(I,xy,uv,tri,ImageSize)


Jr = zeros(ImageSize(1:2));
if(size(I,3)>1), Jg = Jr; Jb = Jr; end

for k=1:size(tri,1)
  
    Vertices=[uv(tri(k,1),:);uv(tri(k,2),:);uv(tri(k,3),:)];
    Vertices2=[xy(tri(k,1),:);xy(tri(k,2),:);xy(tri(k,3),:)];
    
     
    % Get bounding box (ROI)
    boundmin=floor(min(Vertices));
    boundmax=ceil(max(Vertices));
    
    % Bounding box always inside bitmap
    boundmin=max(boundmin,[1 1]);
    boundmax=min(boundmax,[size(Jr,1) size(Jr,2)]);
    
    % The vertices
    p1=Vertices(1,:);
    p2=Vertices(2,:);
    p3=Vertices(3,:);
    
      
    % Normalization factors
    f12 = ( p2(2) - p3(2) ) * p1(1)  + (p3(1) - p2(1) ) * p1(2) + p2(1) * p3(2) - p3(1) *p2(2);
    f20 = ( p3(2) - p1(2) ) * p2(1)  + (p1(1) - p3(1) ) * p2(2) + p3(1) * p1(2) - p1(1) *p3(2);
    f01 = ( p1(2) - p2(2) ) * p3(1)  + (p2(1) - p1(1) ) * p3(2) + p1(1) * p2(2) - p2(1) *p1(2);

    % Lambda Gradient
    g12 = [( p2(2) - p3(2) ) (p3(1) - p2(1) )]/f12;
    g20 = [( p3(2) - p1(2) ) (p1(1) - p3(1) )]/f20;
    g01 = [( p1(2) - p2(2) ) (p2(1) - p1(1) )]/f01;
 
    
    % Center compensation
    c12 = (p2(1) * p3(2) - p3(1) *p2(2))/f12;
    c20 = (p3(1) * p1(2) - p1(1) *p3(2))/f20;
    c01 = (p1(1) * p2(2) - p2(1) *p1(2))/f01;


    [i,j]=ndgrid(boundmin(1):boundmax(1),boundmin(2):boundmax(2));

    % Current location
    r=[i(:) j(:)];
    
    % Interpolation values
    Lambda=[g12(1)*r(:,1)+g12(2)*r(:,2)+c12,  ...
            g20(1)*r(:,1)+g20(2)*r(:,2)+c20,  ...
            g01(1)*r(:,1)+g01(2)*r(:,2)+c01];
    
    % Check if voxel is inside the triangle
    CheckInside=~any((Lambda>1)|(Lambda<0),2);
    
    r(~CheckInside,:)=[];
    Lambda(~CheckInside,:)=[];
    
    posuv=Lambda(:,1)*Vertices2(1,:)+Lambda(:,2)*Vertices2(2,:)+Lambda(:,3)*Vertices2(3,:);
  
       
    ind1=sub2ind([size(Jr,1) size(Jr,2)],r(:,1),r(:,2));
    
    It=image_interpolation_backward(I,posuv-1,'bilinear','replicate',size(posuv(:,1)));
    
    Jr(ind1)=It(:,1);
    if(size(I,3)>1)
        Jg(ind1)=It(:,2);
        Jb(ind1)=It(:,3);
    end
end
J(:,:,1)=Jr;
if(size(I,3)>1)
    J(:,:,2)=Jg;
    J(:,:,3)=Jb;
end

