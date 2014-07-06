function N=ASM_GetContourNormals3D(V,F)
FV.vertices=V;
FV.faces=F;
N=patchnormals(FV);
