function [posV,I_segment]=ASM_ApplyModel3D(Itest,tform,ShapeData,AppearanceData,options)
% Optimalization

% The initial contour is the Mean trainingdata set contour
nl=length(ShapeData.x_mean)/3;
posV=[ShapeData.x_mean(1:nl) ShapeData.x_mean(nl+1:nl*2)  ShapeData.x_mean(nl*2+1:end)];

% Position the initial contour at a location close to the correct location
if(isstruct(tform))
    posV=ASM_align_data_inverse3D(posV,tform);
else
    posV=tform;
end

FV.vertices=posV; FV.faces=ShapeData.Faces;
showcs3(Itest); view(3)
hold on;
h2=[];
handle_fv=patch(FV,'facecolor',[0 0 1],'edgecolor', 'none');   drawnow; pause(1);

% We optimize starting from a rough to a fine image scale
for itt_res=options.nscales:-1:1
    % Scaling of the image
    scale=1 / (2^(itt_res-1));
    
    disp(['scale : ' num2str(scale)]); drawnow;
    PCAData=AppearanceData(itt_res).PCAData;
    
    if((round(scale*1e8)/1e8)==1)
        Itestsmall=Itest;
    else
        %! Warning imgaussian inbouwen!!
        Itestsmall=imresize3d(Itest,scale,[],'cubic','bound');
    end
    
    % Test Starting positions
    if(itt_res==options.nscales&&options.optimizestart)
        PosIn=posV;
        [posVold,Eold]=ASMsearch(Itestsmall,PosIn,ShapeData,AppearanceData,PCAData,options,ShapeData.Faces,itt_res,h2,handle_fv,scale,3);
        
		minp=max(floor(abs(min(PosIn,[],1)-mean(PosIn,1)))+1,floor(mean(PosIn,1)-mean(size(Itest))*0.15));
		maxp=min(size(Itest)-floor(abs(max(PosIn,[],1)-mean(PosIn,1))),ceil(mean(PosIn,1)+mean(size(Itest))*0.15));
		[x,y,z]=ndgrid(linspace(minp(1),maxp(1),3),linspace(minp(2),maxp(2),3),linspace(minp(3),maxp(3),3));
		offset=[x(:) y(:) z(:)];
        offset=bsxfun(@minus,offset,mean(PosIn,1));
     
        verbose=options.verbose;
        options.verbose=false;
		disp(['Number of start positions: ' num2str(size(offset,1))]); drawnow;
        for i=1:size(offset,1);
            disp(['start positions: ' num2str(i)]); drawnow;
  
            posV=PosIn;
            posV(:,1)=posV(:,1)+offset(i,1);
            posV(:,2)=posV(:,2)+offset(i,2);
            posV(:,3)=posV(:,3)+offset(i,3);
            [posV,E]=ASMsearch(Itestsmall,posV,ShapeData,AppearanceData,PCAData,options,ShapeData.Faces,itt_res,h2,handle_fv,scale,3);
            if(E<Eold), Eold=E; posVold=posV; end
        end
        options.verbose=verbose;
        posV=posVold;
    end
    
    % Do 50 ASM itterations
    nsearch=options.nsearch(min(itt_res,length(options.nsearch)));
    [posV,E]=ASMsearch(Itestsmall,posV,ShapeData,AppearanceData,PCAData,options,ShapeData.Faces,itt_res,h2,handle_fv,scale,nsearch);
    disp(['Current Error : ' num2str(E)]);
end
if(nargout>1)
    if(isempty(ShapeData.Faces))
        I_segment = [];
    else
        I_segment = polygon2voxel(struct('vertices',posV,'Faces',ShapeData.Faces),size(Itest),'clamp', false);
        I_segment = imfill(I_segment,'holes');
    end
end




function f=calculateMovementEnergyImage(gt,dgt,k2,ns2,nl,originalsearch,AppearanceData,PCAData,itt_res)
% Loop through all contour points
f=zeros(ns2,nl);
[x,y]=ndgrid(1:k2,1:ns2);  indGI=x+y-1;

for j=1:nl
    drawnow
    % Search through the large sampeled profile, for the optimal
    % position
    % A profile from the large profile, with the length
    % of the trainingprofile (for rgb image 3x as long)
    ind=sub2ind(size(gt),indGI(:),repmat(j,numel(indGI),1));
    if(originalsearch)
        gi=reshape(dgt(ind),k2,ns2);
        % Calculate the Mahalanobis distance from the current
        % profile, to the training data sets profiles through
        % an inverse correlation matrix.
        for i=1:size(gi,2)
            v=(gi(:,i)- AppearanceData(itt_res).Landmarks(j).dg_mean);
            f(i,j)=v'*AppearanceData(itt_res).Landmarks(j).Sinv*v;
        end
    else
        gi=reshape(gt(ind),k2,ns2);
        % Calculate the PCA parameters, and normalize them
        % with the variances.
        % (Seems to work better with color images than
        % the original method)
        bc = PCAData(j).Evectors'*(gi-repmat(PCAData(j).Emean,1,ns2));
        bc = bc./repmat(sqrt(PCAData(j).Evalues),1,ns2);
        f(:,j)=sum(bc.^2,1);
    end
    % Find the lowest Mahalanobis distance, and store it
    % as movement step
end

function [posV,E]=ASMsearch(Itestsmall,posV,ShapeData,AppearanceData,PCAData,options,Faces,itt_res,h2,handle_fv,scale,nsearch)
nl=length(ShapeData.x_mean)/3;
for itt=1:nsearch
    FVR.vertices=posV; FVR.faces=Faces;
    disp(['itteration : ' num2str(itt)]); drawnow;
    
    % Calculate the normals of the contour points.
    N=ASM_GetContourNormals3D(posV,Faces);
    
    % Create long intensity profiles on the contour normals, for search
    % of best point fit, using the correlation matrices made in the
    % appearance model
    
    % Total Intensity line-profile needed, is the traininglength + the
    % search length in both normal directions
    n = options.k + options.ns;
    
    % Get the intensity profiles of all landmarks and there first order
    % derivatives
    [gt,dgt]=ASM_getProfileAndDerivatives3D(Itestsmall,posV*scale,N,n);
    
    % Calculate optimal movement
    k2=2*options.k+1;
    ns2=options.ns*2+1;
    originalsearch=options.originalsearch;
    f=calculateMovementEnergyImage(gt,dgt,k2,ns2,nl,originalsearch,AppearanceData,PCAData,itt_res);
    
    % Calculate regularization error
    E=sum(f(options.ns,:));
	if(options.verbose)
		disp(['Cf. Error : ' num2str(E)]);
    end
	
    [temp,i]=min(f);
    movement=((i-1)-options.ns);
    
    % Move the points to there new optimal positions
    posV=posV+(1/scale)*repmat(movement',1,3).*N;
    
    FVR.vertices=posV; FVR.faces=Faces;
    
    % Show the new positions
    if(options.verbose)
		if(ishandle(h2)), delete(h2); end
		h2=plot3(posV(:,1),posV(:,2),posV(:,3),'r.'); drawnow('expose');
	end
	
    % Outlier removal
    posV=OutlierCorrection(posV,ShapeData,Faces);
    
    % Remove translation and rotation, as done when training the
    % model.
    
    posVIn=posV;
    [posV,tform]=ASM_align_data3D(posVIn,ShapeData.MeanVertices);
    
    % Describe the model by a vector b with model parameters
    x_search=[posV(:,1);posV(:,2);posV(:,3)];
    b = ShapeData.Evectors'*(x_search-ShapeData.x_mean);
    
    % Limit the model parameters based on what is considered a nomal
    % contour, using the eigenvalues of the PCA-model
    maxb=options.m*sqrt(ShapeData.Evalues);
    b=max(min(b,maxb),-maxb);
    
    % Transform the model parameter vector b, back to contour positions
    x_search = ShapeData.x_mean + ShapeData.Evectors*b;
    posV(:,1)=x_search(1:nl)';
    posV(:,2)=x_search(nl+1:nl*2)';
    posV(:,3)=x_search(nl*2+1:end)';
    
    % Now add the previously removed translation and rotation
    posV=ASM_align_data_inverse3D(posV,tform);
    
    % Calculate regularization error
    %E=sum(sqrt(sum((posV-posVIn).^2,2)));
    %disp(['C. Error : ' num2str(E)]);
    
    if(itt_res==1)
        %   [~,~,posV]=point_registration(size(Itest),posV,posVIn,struct('Verbose',false,'MaxRef',2));
    end
    
    if(options.verbose)
        if(ishandle(handle_fv)), delete(handle_fv); end
        FV.vertices=posV; FV.faces=Faces;
        handle_fv=patch(FV,'facecolor',[0 1 0],'edgecolor', 'none'); drawnow; pause(1);
    end
end

function posIn=OutlierCorrection(pos,ShapeData,Faces)
posIn=pos;
pos=ASM_align_data3D(posIn,ShapeData.MeanVertices);
x_search=[pos(:,1);pos(:,2);pos(:,3)];
V=abs(bsxfun(@times,ShapeData.Evectors,x_search-ShapeData.x_mean));
V=reshape(V,size(pos,1),[]);
V=max(bsxfun(@rdivide,V,sum(V,1)),[],2);
Outliers=V>mean(V,1)*2.5;
index=find(Outliers);
if(~isempty(index))
    % Set Outliers to mean of neighbors
    for i=1:length(index)
        F=Faces(any(Faces==index(i),2),:);
        F=unique(F(:));
        F(F==index(i))=[];
        V(index(i),:)=mean(V(F,:),1);
    end
end