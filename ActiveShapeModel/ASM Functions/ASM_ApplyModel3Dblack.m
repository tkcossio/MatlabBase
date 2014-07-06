function posV=ASM_ApplyModel3Dblack(Itest,tform,ShapeData,AppearanceData,options,Faces)
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

FV.vertices=posV; FV.faces=Faces;
showcs3(Itest); view(3)
hold on;
h2=[];
handle_fv=patch(FV,'facecolor',[0 0 1],'edgecolor', 'none');   drawnow; pause(1);

pd=0;
% We optimize starting from a rough to a fine image scale
for itt_res=options.nscales:-1:1
    % Scaling of the image
    scale=1 / (2^(itt_res-1));

    disp(['scale : ' num2str(scale)]); drawnow;
    
    PCAData=AppearanceData(itt_res).PCAData;
     
    % Do 50 ASM itterations
    if((round(scale*1e8)/1e8)==1)
        Itestsmall=Itest;
    else
    	Itestsmall=imresize3d(Itest,scale,[],'cubic','bound');
    end
    for itt=1:options.nsearch(min(itt_res,length(options.nsearch)))
        pd=pd+1; 
        FVR.vertices=posV; FVR.faces=Faces;
        %save(['demoasm' num2str(pd)],'FVR');
        disp(['itteration : ' num2str(itt)]); drawnow;
           
        % Calculate the normals of the contour points.
        N=ASM_GetContourNormals3D(posV,Faces);

        % Create long intensity profiles on the contour normals, for search
        % of best point fit, using the correlation matrices made in the
        % appearance model
        
        % Total Intensity line-profile needed, is the search length 
        % in both normal directions
        n =  options.ns;

        % Get the intensity profiles of all landmarks and there first order
        % derivatives
        gt=ASM_getProfileAndDerivatives3D(Itestsmall,(posV-0.5)*scale+0.5,N,n);

		% Calculate optimal movement
		 
        [temp,i]=min(imfilter(gt,ones(4,1)));
        movement=((i-1)-options.ns);

        % Move the points to there new optimal positions
        posV=posV+(1/scale)*repmat(movement',1,3).*N;
        
        FVR.vertices=posV; FVR.faces=Faces;
        %save(['demoasms' num2str(pd)],'FVR');
        
        % Show the new positions
        if(ishandle(h2)), delete(h2); end
        h2=plot3(posV(:,1),posV(:,2),posV(:,3),'r.'); drawnow('expose'); 

        % Remove translation and rotation, as done when training the
        % model.
        [posV,tform]=ASM_align_data3D(posV,ShapeData.MeanVertices);
        
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
        if(ishandle(handle_fv)), delete(handle_fv); end
		FV.vertices=posV; FV.faces=Faces;
        handle_fv=patch(FV,'facecolor',[0 1 0],'edgecolor', 'none'); drawnow; pause(1);
    end
end





