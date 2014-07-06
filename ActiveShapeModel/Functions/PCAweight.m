function [Evalues, Evectors1, Evectors2, x_mean]=PCAweight(x,w)
% Weighted PCA using Single Value Decomposition
% Obtaining mean vector, eigenvectors and eigenvalues
%
% [Evalues, Evectors1,Evectors2, x_mean]=PCA(x,w);
%
% inputs,
%   x : M x N matrix with M the trainingvector length and N the number
%              of training data sets
%   w : M x 1 weights of the features
%
% outputs,
%   Evalues : The eigen values of the data
%   Evector1 : The forward, eigen vectors of the data
%   Evector2 : The backward, eigen vectors of the data
%   x_mean : The mean training vector
%
% Example,
%  % TrainingData
%  x=rand(10,6);
%  w=[5 1 1 1 1 1 1 1 1 1]';
% 
%  % Do the weighted PCA
%  [Evalues, Evectors1, Evectors2, x_mean]=PCAweight(x,w);
%
%  % Testdata
%  y=rand(10,1);
%  
%  % From normal values to PCA model values
%  b = Evectors1*(y-x_mean);
%
%  % Limit the values 
%  maxb=3*sqrt(Evalues);
%  b=max(min(b,maxb),-maxb);
%
%  % Go back to normal values
%  ym=x_mean + Evectors2*b;
%
%  disp(y-ym)
%
%
%

% Use the weights to change the variance of some features
x=bsxfun(@times,x,w);

s=size(x,2);
% Calculate the mean 
x_mean=sum(x,2)/s;

% Substract the mean
x2=(x-repmat(x_mean,1,s))/ sqrt(s-1);

% Do the SVD 
[U2,S2] = svd(x2,0);

% Calculate the eigenvectors and values of the weighted data
Evalues=diag(S2).^2;
Evectors=bsxfun(@times,U2,sign(U2(1,:)));

% Remove weighting from mean
winv=ones(size(w))./w;
x_mean=x_mean.*winv;

% Forward and Backward Evectors including the weights
Evectors1=bsxfun(@times,Evectors,w)';
Evectors2=bsxfun(@times,Evectors,winv);
