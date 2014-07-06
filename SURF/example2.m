% Example 2, Corresponding points
% Load images
  I1 = phantom(512,512);
I2 = phantom(256,256);
I2(257:512,1:512) = 0;
% Get the Key Points
  Options.upright=true;
  Options.tresh=0.0001;
  Ipts1=OpenSurf(I1,Options);
  Ipts2=OpenSurf(I2,Options);
% Put the landmark descriptors in a matrix
  D1 = reshape([Ipts1.descriptor],64,[]); 
  D2 = reshape([Ipts2.descriptor],64,[]); 
% Find the best matches
  err=zeros(1,length(Ipts1));
  cor1=1:length(Ipts1); 
  cor2=zeros(1,length(Ipts1));
  for i=1:length(Ipts1),
      distance=sum((D2-repmat(D1(:,i),[1 length(Ipts2)])).^2,1);
      [err(i),cor2(i)]=min(distance);
  end
% Sort matches on vector distance
  [err, ind]=sort(err); 
  cor1=cor1(ind); 
  cor2=cor2(ind);
% Show both images
  I = zeros([size(I1,1) size(I1,2)*2 size(I1,3)]);
  I(:,1:size(I1,2),:)=I1; I(:,size(I1,2)+1:size(I1,2)+size(I2,2),:)=I2;
  figure, imagesc(I); hold on;
% Show the best matches
  for i=1:30,
      c=rand(1,3);
      plot([Ipts1(cor1(i)).x Ipts2(cor2(i)).x+size(I1,2)],[Ipts1(cor1(i)).y Ipts2(cor2(i)).y],'-','Color',c)
      plot([Ipts1(cor1(i)).x Ipts2(cor2(i)).x+size(I1,2)],[Ipts1(cor1(i)).y Ipts2(cor2(i)).y],'o','Color',c)
  end
  
%%
%Define base image
I1 = imread('C:\Users\cossitk1\Documents\ToyData\Kimia216\bird02.pgm');
I2 = imread('C:\Users\cossitk1\Documents\ToyData\Kimia216\bird04.pgm');
MAX_ROWS = max([size(I1,1), size(I2,1)]);
MAX_COLS = max([size(I2,2), size(I2,2)]);

I1(end:MAX_ROWS, :) = 0;
I1(:, end:MAX_COLS) = 0;
I2(end:MAX_ROWS, :) = 0;
I2(:, end:MAX_COLS) = 0;

I1 = double(I1);
I2 = double(I2);

%I1 = phantom(512,512);
%I2 = phantom(256,256);
% Get the Key Points
  Options.upright=true;
  Options.tresh=0.0001;
  Ipts1=OpenSurf(I1,Options);
  Ipts2=OpenSurf(I2,Options);
% Put the landmark descriptors in a matrix
  D1 = reshape([Ipts1.descriptor],64,[]); 
  D2 = reshape([Ipts2.descriptor],64,[]); 
% Find the best matches
  err=zeros(1,length(Ipts1));
  cor1=1:length(Ipts1); 
  cor2=zeros(1,length(Ipts1));
  for i=1:length(Ipts1),
      distance=sum((D2-repmat(D1(:,i),[1 length(Ipts2)])).^2,1);
      [err(i),cor2(i)]=min(distance);
  end
% Sort matches on vector distance
  [err, ind]=sort(err); 
  cor1=cor1(ind); 
  cor2=cor2(ind);
% Show both images
  I = zeros([size(I1,1) size(I1,2)*2 size(I1,3)]);
  I(:,1:size(I1,2),:)=I1; I(:,size(I1,2)+1:size(I1,2)+size(I2,2),:)=I2;
  figure, imagesc(I); hold on;
% Show the best matches

err_match = find(err < 0.2);

  for match_idx=1:length(err_match),
      i = err_match(match_idx);
      c=rand(1,3);
      plot([Ipts1(cor1(i)).x Ipts2(cor2(i)).x+size(I1,2)],[Ipts1(cor1(i)).y Ipts2(cor2(i)).y],'-','Color',c)
      plot([Ipts1(cor1(i)).x Ipts2(cor2(i)).x+size(I1,2)],[Ipts1(cor1(i)).y Ipts2(cor2(i)).y],'o','Color',c)
  end
 
 