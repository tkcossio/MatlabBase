function N=ASM_GetContourNormals2D(V,L)
% This function calculates the normals, of the contour points
% using the neighbouring points of each contour point, and 
% forward an backward differences on the end points
%
% N=GetContourNormals(V,L)
%
% inputs,
%   V : List of Vertices 2 x N
%   L : Line list, with indices to the vertices 2 x M
%
% outputs,
%   N : The normals of the Vertices
%

% Derivatives of contour
DT=V(L(:,1),:)-V(L(:,2),:);
D1=zeros(size(V));
D2=zeros(size(V));
D1(L(:,1),:)=DT;
D2(L(:,2),:)=DT;
D=D1+D2;
L=sqrt(D(:,1).^2+D(:,2).^2);
N(:,1)= D(:,2)./L;
N(:,2)=-D(:,1)./L;
