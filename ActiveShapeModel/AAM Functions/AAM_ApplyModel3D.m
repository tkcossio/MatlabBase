function [pos,I_model,I_segment]=AAM_ApplyModel3D(ItestLarge,tformLarge,Data,options)
pos = InitialShape(tformLarge,Data,options);

% Loop through the 4 image size scales
for scale=options.nscales:-1:1
    % Get the PCA model for this scale
    R=Data{scale}.R;
    ShapeAppearanceData=Data{scale}.ShapeAppearanceData;
    ShapeData=Data{scale}.ShapeData;
    AppearanceData=Data{scale}.AppearanceData;
    
    % The image scaling of the scale-itteration
    scaling=2^(-(scale-1));
    pos=(pos)*scaling;
        
    % Transform the image and coordinates offset, to the current scale
    if(round(scaling*1000)~=1000)
        Itest=imresize3d(ItestLarge,scaling,[],'cubic','replicate');
    else
        Itest=ItestLarge;
    end
    
    % Start a new figure
    % show current test image, and initial contour
    % From real image coordinates to -> algined coordinates
    
    FV.vertices=pos;
    FV.faces=ShapeData.Faces;
    if(options.verbose)
        showcs3(Itest); view(3)
        hold on;
        if(isempty(FV.faces))
            plot3(FV.vertices(:,1),FV.vertices(:,2)  ,FV.vertices(:,3) );  drawnow; pause(1);
        else
            patch(FV,'facecolor',[0 0 1],'edgecolor', 'none');   drawnow; pause(1);
        end
    end
    
    % Convert the inital surface shape to Model and Pose parameters
    cc=InitShapeToModelPoseParameters(pos,Itest,options,ShapeData,AppearanceData,ShapeAppearanceData);
    
    % Try several init positions
    if(scale==options.nscales&&options.optimizestart)
        PosIn=pos;
        optionsIn=options; options.nsearch=3;
        Eold=ErrorModelImage(cc,AppearanceData,ShapeAppearanceData,ShapeData,Itest,options);
        ccold=cc;
 
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
            pos=PosIn;
            pos(:,1)=pos(:,1)+offset(i,1);
            pos(:,2)=pos(:,2)+offset(i,2);
            pos(:,3)=pos(:,3)+offset(i,3);
            [cc,E]=SearchAAM(InitShapeToModelPoseParameters(pos,Itest,options,ShapeData,AppearanceData,ShapeAppearanceData),Itest,options,ShapeData,AppearanceData,ShapeAppearanceData,R);
            if(E<Eold),  Eold=E; ccold=cc; end
        end
        cc=ccold; options=optionsIn;    
        options.verbose=verbose;
            % Adjusted distance
        %disp(['Offset Adjust distance in pixels' num2str(dis)]);
    end
    
    % Starting step size
    if(options.optimizecmiddle)
        [cc,E,pos]=SearchAAMsimplex(cc,Itest,options,ShapeData,AppearanceData,ShapeAppearanceData);
    else
        [cc,E,pos]=SearchAAM(cc,Itest,options,ShapeData,AppearanceData,ShapeAppearanceData,R);
    end
    
    hold on;
    FV.vertices=pos; FV.faces=ShapeData.Faces;
    if(options.verbose)
        if(isempty(FV.faces))
            plot3(FV.vertices(:,1),FV.vertices(:,2)  ,FV.vertices(:,3) ); drawnow; pause(1);
        else
            patch(FV,'facecolor',[1 0 0],'edgecolor', 'none','facealpha',0.5); drawnow; pause(1);
        end
    end
    
    disp(['Current Error : ' num2str(E)])
    
    if(options.optimizecend)
        cc = fminsearch(@(x)ErrorModelImage(x,AppearanceData,ShapeAppearanceData,ShapeData,Itest,options),cc,struct('MaxIter',2*options.nsearch*ceil(2^(scale-1)),'Display','off'));
        [E,~,~,pos]=ErrorModelImage(cc,AppearanceData,ShapeAppearanceData,ShapeData,Itest,options);
        disp(['Current Error : ' num2str(E)])
    end
    
    hold on;
    FV.vertices=pos; FV.faces=ShapeData.Faces;
    if(options.verbose)
        if(isempty(FV.faces))
            plot3(FV.vertices(:,1),FV.vertices(:,2)  ,FV.vertices(:,3) ); drawnow; pause(1);
        else
            patch(FV,'facecolor',[1 0 1],'edgecolor', 'none','facealpha',0.5); drawnow; pause(1);
        end
    end
    
    % Scale the contour for the next itteration
    pos=(pos)/scaling;
end
[~,g,g_model,pos]=ErrorModelImage(cc,AppearanceData,ShapeAppearanceData,ShapeData,Itest,options);

% Convert the model vector to model-image-volume
c1=std(g(:)); c2=mean(g(:));
I_texture=AAM_Vector2Appearance3D(g_model*c1+c2,AppearanceData.ObjectPixels,ShapeData.TextureSize);

% Warp the model-image-volume to right location in the the testvolume
input_points= AppearanceData.base_points;
base_points=pos;
xyz=[input_points(:,2) input_points(:,1) input_points(:,3)];
uvw=[base_points(:,2) base_points(:,1) base_points(:,3)];
I_model = warp_tetrahedron(I_texture,xyz,uvw, size(ItestLarge));

% Make the segmentation
if(isempty(ShapeData.Faces))
    I_segment = [];
else
    I_segment = polygon2voxel(struct('vertices',pos,'faces',ShapeData.Faces),size(ItestLarge),'none', false);
    I_segment = imfill(I_segment,'holes');
end


function [cc,E,pos]=SearchAAMsimplex(cc,Itest,options,ShapeData,AppearanceData,ShapeAppearanceData)
for k=1:2
    cc = fminsearch(@(x)ErrorModelImage(x,AppearanceData,ShapeAppearanceData,ShapeData,Itest,options),cc,struct('MaxIter',8*options.nsearch*ceil(2^(scale-1)),'Display','off'));
    [E,~,~,pos]=ErrorModelImage(cc,AppearanceData,ShapeAppearanceData,ShapeData,Itest,options);
end

function [cc,E,pos]=SearchAAM(cc,Itest,options,ShapeData,AppearanceData,ShapeAppearanceData,R)
if(options.scale3), posen=10; else posen=7; end
w=1;
cc_old=[]; pos_old=[]; Eold=inf;
% Search Itterations
for i=1:options.nsearch
    [E,g,g_model,pos]=ErrorModelImage(cc,AppearanceData,ShapeAppearanceData,ShapeData,Itest,options);
    
    % Go back to the old location of the previous itteration, if the
    % error was lower.
    if(E>Eold)
        % Not always go back if the error becomes higher, sometimes
        % stay on the higher error (like in simulated annealing)
        % Try a smaller stepsize
        w=w*0.9;
        cc=cc_old; pos=pos_old;
    else
        w=w*1.1;  Eold=E;
    end
    
    % Store model /pose parameters for next itteration
    cc_old=cc;
    pos_old=pos;
    
    % Calculate the needed model parameter update using the
    % search model matrix
    cc_diff=R*(g-g_model);
    
    % Update the ShapeApppearance Parameters,
    % and stay within 3 (m) standard deviations
    cc(1:end-posen)=cc(1:end-posen)+cc_diff(1:end-posen)*w;
    cc(1:end-posen)=max(min(cc(1:end-posen),options.m),-options.m);
    
    % Update the pose parameters
    cc(end-posen+1:end)=cc(end-posen+1:end)+cc_diff(end-posen+1:end)*w;
end
if(E>Eold), cc=cc_old;  pos=pos_old; E=Eold; end



function cc=InitShapeToModelPoseParameters(pos,Itest,options,ShapeData,AppearanceData,ShapeAppearanceData)
[pos_align, tform]=AAM_align_data3D(pos, ShapeData.MeanVertices,options);
x=[pos_align(:,1);pos_align(:,2);pos_align(:,3)];

% Sample the image intensities
g=AAM_Appearance2Vector3D(Itest,pos, AppearanceData.base_points, AppearanceData.ObjectPixels,ShapeData.TextureSize,AppearanceData.Tetra,options,ShapeData.Faces);
g=AAM_NormalizeAppearance3D(g,options);

% Go from image intesities and contour to ShapeAppearance parameters
b1 = ShapeAppearanceData.Ws * ShapeData.Evectors' * (x-ShapeData.x_mean);
b2 = AppearanceData.Evectors' * (g-AppearanceData.g_mean); b2=double(b2);
b = [b1;b2];
cin = ShapeAppearanceData.Evectors'*(b -ShapeAppearanceData.b_mean);

x2=x;
maxc=options.m*sqrt(ShapeAppearanceData.Evalues);
if(options.optimizecstart)
    c = lsqnonlin(@(x)costc(x,double(x2),ShapeAppearanceData,ShapeData),double(cin),-maxc,maxc,optimset('Display','off','MaxIter',10));
else
    c = cin;
end

%  Storage of ShapeAppeanance parameters, pose parameters, last
%  error between model and intensities. Used if old location
%  had a smaller intensity error.

cc=CombineModelPoseParameters(tform,options,c,ShapeData,ShapeAppearanceData);


function pos = InitialShape(tformLarge,Data,options)
% We start at the coarse scale
scale=options.nscales; scaling=2^(-(scale-1));

% Get the PCA model for this scale
ShapeAppearanceData=Data{scale}.ShapeAppearanceData;
ShapeData=Data{scale}.ShapeData;

% Use the mean ShapeAppearance parameters to go get an initial contour
b = ShapeAppearanceData.b_mean;
b1 = b(1:(length(ShapeAppearanceData.Ws)));
b1= ShapeAppearanceData.Ws\b1;

if(isstruct(tformLarge))
    % Initial (mean) aligned coordinates
    x = ShapeData.x_mean + ShapeData.Evectors*b1;
    
    % The real image coordinates
    nl=length(x)/3;
    pos(:,1)=x(1:nl);
    pos(:,2)=x(nl+1:nl*2);
    pos(:,3)=x(nl*2+1:end);
    
    % Transform the coordinates to match the coarse scale
    tform=tformLarge;
    %tform.offsetv=(tform.offsetv-0.5)*scaling+0.5;
    tform.offsetv=tform.offsetv*scaling;
    pos=AAM_align_data_inverse3D(pos,tform,options);
    pos=(pos)/scaling;
else
    pos=tformLarge;
end

function E=costc(c,x2,ShapeAppearanceData,ShapeData)
% Probably also solve able with lsqlin (linear bound solve)...

% Go from ShapeAppearance Parameters to aligned shape coordinates
b = ShapeAppearanceData.b_mean + ShapeAppearanceData.Evectors*c;
b1 = b(1:(length(ShapeAppearanceData.Ws)));
b1= ShapeAppearanceData.Ws\b1;
x = ShapeData.x_mean + ShapeData.Evectors*b1;
E=double(x(:)-x2(:));


function [E,g,g_model,pos]=ErrorModelImage(cc,AppearanceData,ShapeAppearanceData,ShapeData,Itest,options)
% Split Model + Pose parameters
[tform,c]=SplitModelPoseParameters(options,cc,ShapeData,ShapeAppearanceData);

% Go from ShapeAppearance Parameters to aligned shape coordinates
b = ShapeAppearanceData.b_mean + ShapeAppearanceData.Evectors*c;
b1 = b(1:(length(ShapeAppearanceData.Ws)));
b1= ShapeAppearanceData.Ws\b1;
x = ShapeData.x_mean + ShapeData.Evectors*b1;

% From aligned coordinates to real image coordinates
nl=length(x)/3;
pos=AAM_align_data_inverse3D([x(1:nl) x(nl+1:nl*2) x(nl*2+1:end)],tform,options);

% Sample the intensities
g=AAM_Appearance2Vector3D(Itest,pos, AppearanceData.base_points, AppearanceData.ObjectPixels,ShapeData.TextureSize,AppearanceData.Tetra,options,ShapeData.Faces);
g=AAM_NormalizeAppearance3D(g,options);

% Go from intensities and shape back to ShapeAppearance Parameters
b1 = ShapeAppearanceData.Ws * ShapeData.Evectors' * (x-ShapeData.x_mean);
b2 = AppearanceData.Evectors' * (g-AppearanceData.g_mean); b2=double(b2);
b = [b1;b2];
c2 = ShapeAppearanceData.Evectors'*(b -ShapeAppearanceData.b_mean);
if(options.usemodelc), c2=c;end

% Go from ShapeAppearance Parameters back to model intensities
b = ShapeAppearanceData.b_mean + ShapeAppearanceData.Evectors*c2;
b2 = b(size(ShapeAppearanceData.Ws,1)+1:end);
g_model = AppearanceData.g_mean + AppearanceData.Evectors*b2;

% Difference between model and real image intensities
if(options.logerror)
    E=sum(log(1+(g-g_model).^2/2*options.sigmas^2));
else
    E=sum((g-g_model).^2);
end
if(isnan(E)), E=1e100; end


function cc=CombineModelPoseParameters(tform,options,c,ShapeData,ShapeAppearanceData)
if(options.scale3), posen=10; else posen=7; end
cc=[c(:);zeros(posen,1)];
cc(1:end-posen)=cc(1:end-posen)./sqrt(ShapeAppearanceData.Evalues);

if(~options.scale3),
    cc(end-6:end-4)=tform.offsetv;
    cc(end-3:end)=tform.offsetq;
else
    cc(end-9:end-7)=tform.offsetv;
    cc(end-6:end-3)=tform.offsetq;
    cc(end-2:end)=tform.offsets;
end

if(options.posevariance)
    if(~options.scale3)
        cc(end-posen+1:end)=cc(end-posen+1:end)./[ShapeData.TVariance(:);ShapeData.QVariance(:)];
    else
        cc(end-posen+1:end)=cc(end-posen+1:end)./[ShapeData.TVariance(:);ShapeData.QVariance(:);ShapeData.SVariance(:)];
    end
end

function [tform,c]=SplitModelPoseParameters(options,cc,ShapeData,ShapeAppearanceData)
if(options.scale3), posen=10; else posen=7; end

% Scale Pose and Model parameters back (with variance)
c=cc(1:end-posen).*sqrt(ShapeAppearanceData.Evalues);
if(options.posevariance)
    if(~options.scale3)
        cc(end-posen+1:end)=cc(end-posen+1:end).*[ShapeData.TVariance(:);ShapeData.QVariance(:)];
    else
        cc(end-posen+1:end)=cc(end-posen+1:end).*[ShapeData.TVariance(:);ShapeData.QVariance(:);ShapeData.SVariance(:)];
    end
end
cc(1:end-posen)=max(min(cc(1:end-posen),options.m),-options.m);

if(~options.scale3),
    tform.offsetv=cc(end-6:end-4);
    tform.offsetq=cc(end-3:end);
else
    tform.offsetv=cc(end-9:end-7);
    tform.offsetq=cc(end-6:end-3);
    tform.offsets=cc(end-2:end);
end

