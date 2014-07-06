function J=warp_tetrahedron_double(I,xyz,uvw,tetra,ImageSize)


J = zeros(ImageSize(1:3));

for q=1:size(tetra,1)
    Vertices =[uvw(tetra(q,1),:);uvw(tetra(q,2),:);uvw(tetra(q,3),:);uvw(tetra(q,4),:)];
    Vertices2=[xyz(tetra(q,1),:);xyz(tetra(q,2),:);xyz(tetra(q,3),:);xyz(tetra(q,4),:)];
    
    % Get bounding box (ROI)
    boundmin=floor(min(Vertices));
    boundmax=ceil(max(Vertices));
    
    % Bounding box always inside bitmap
    boundmin=max(boundmin,[1 1 1]);
    boundmax=min(boundmax,[size(J,1) size(J,2) size(J,3)]);
    
    % The vertices
    for kq=1:1
        switch(kq)
            case 1
                ind=[1 2 3 4];
            case 2
                ind=[4 1 2 3];
            case 3
                ind=[3 4 1 2];
            case 4
                ind=[2 3 4 1];
        end
                
        p1=Vertices(ind(1),:);
        p2=Vertices(ind(2),:);
        p3=Vertices(ind(3),:);
        p4=Vertices(ind(4),:);

        % Define matrix T
        T=[p1(:)-p4(:), p2(:)-p4(:),p3(:)-p4(:)];
        detT=T(1,1)*(T(2,2)*T(3,3)-T(2,3)*T(3,2))+ ...
             T(1,2)*(T(2,3)*T(3,1)-T(3,3)*T(2,1))+ ...
             T(1,3)*(T(2,1)*T(3,2)-T(2,2)*T(3,1));

     
        Tinv=(1/detT)*[cross(T(:,2),T(:,3)),cross(T(:,3),T(:,1)),cross(T(:,1),T(:,2))]';
           
        % Center compensation
        c1=-Tinv(1,:)*p4(:);
        c2=-Tinv(2,:)*p4(:);
        c3=-Tinv(3,:)*p4(:);
    
        [i,j,k]=ndgrid(boundmin(1):boundmax(1),boundmin(2):boundmax(2),boundmin(3):boundmax(3));

        % Current location
        r=[i(:) j(:) k(:)];

        % Interpolation values
        Lambda=[Tinv(1,1)*r(:,1)+Tinv(1,2)*r(:,2)+Tinv(1,3)*r(:,3)+c1, ...
                Tinv(2,1)*r(:,1)+Tinv(2,2)*r(:,2)+Tinv(2,3)*r(:,3)+c2, ...
                Tinv(3,1)*r(:,1)+Tinv(3,2)*r(:,2)+Tinv(3,3)*r(:,3)+c3];
        Lambda(:,4)=ones(size(Lambda(:,1)))-Lambda(:,1)-Lambda(:,2)-Lambda(:,3);
        
            switch(kq)
            case 1
                Lambdat=Lambda;
            case 2
                Lambdat=Lambdat+Lambda(:,[2 3 4 1]);
            case 3
                Lambdat=Lambdat+Lambda(:,[3 4 1 2]);
            case 4
                Lambdat=Lambdat+Lambda(:,[4 1 2 3]);
            end
    end
    Lambda=Lambdat/kq;
       
    % Check if voxel is inside tetrahedron
    CheckInside=~any((Lambda>1.0001)|(Lambda<-0.0001),2);
    r(~CheckInside,:)=[];
    Lambda(~CheckInside,:)=[];
    
    posuvw=Lambda(:,1)*Vertices2(1,:)+Lambda(:,2)*Vertices2(2,:)+Lambda(:,3)*Vertices2(3,:)+Lambda(:,4)*Vertices2(4,:);
    
      
    ind1=sub2ind([size(J,1) size(J,2) size(J,3)],r(:,1),r(:,2),r(:,3));
    
    It=image_interpolation_backward(I,posuvw-1,'bilinear','replicate',size(posuvw(:,1)));
    J(ind1)=It(:);
end
