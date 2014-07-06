#include "mex.h"
#include "math.h"
#include "image_interpolation.h"
#include "multiple_os_thread.h"
#include "string.h"

/* The matlab mex function */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[] )
{
    double *Iout, *Iin, *XY, *UV, *TRI, *SizeO;
    
    const mwSize *Idims;
    const mwSize *Tdims;
    const mwSize *Vdims;
    
    int Isize[3] = {1,1,1};
    int Osize[3] = {1,1,1};
    int Tsize[2] = {1,1};
    int Vsize[2] = {1,1};
    int Oslice;
    /* loop variable */
    int q, i, j, rgb;
    
    /* Current Vertices indices*/
    int t0, t1, t2;
    
    /* Current voxel/pixel */
    double Ipixel[3]={0,0,0};
    
    /* Bounding box polygon */
    int boundmin[2];
    int boundmax[2];
    
    /* Vertices */
    double p0[2],p1[2], p2[2];
    double x0[2],x1[2], x2[2];
    
    /* Barycentric variables */
    double f12, f20, f01;
    double g12[2],g20[2],g01[2];
    double c12,c20,c01;
    double posuv[2];
    
    /* Interpolation percentages */
    double Lambda[3];
    double Lambdat[3];
    
    /* Boundary values */
    const double mval=-0.0001;
    const double bval=1.0001;
    
    
    
    /* function J=warp_triangle_double(I,xy,uv,tri,ImageSize) */
    Iin=mxGetPr(prhs[0]);
    XY=mxGetPr(prhs[1]);
    UV=mxGetPr(prhs[2]);
    TRI=mxGetPr(prhs[3]);
    SizeO=mxGetPr(prhs[4]);
    
    /* Input Image size */
    Idims = mxGetDimensions(prhs[0]);
    Isize[0] = Idims[0];
    Isize[1] = Idims[1];
    
    /* Input Number of Polygons */
    Tdims = mxGetDimensions(prhs[3]);
    Tsize[0] = Tdims[0];
    Tsize[1] = Tdims[1];
    
    /* Input Number of Polygons */
    Vdims = mxGetDimensions(prhs[1]);
    Vsize[0] = Vdims[0];
    Vsize[1] = Vdims[1];
    
    
    if(mxGetNumberOfDimensions(prhs[0])>2) { Isize[2] = 3; }
    
    /* Output Image size */
    Osize[0]=(int)SizeO[0];
    Osize[1]=(int)SizeO[1];
    Osize[2]=Isize[2];
    Oslice=Osize[0]*Osize[1];
    
    /* Create empty array for output */
    plhs[0] = mxCreateNumericArray(3, Osize, mxDOUBLE_CLASS, mxREAL);
    Iout=mxGetPr(plhs[0]);
    if(nrhs<6)
    {
        Iout=memset(Iout,0,Osize[0]*Osize[1]*Osize[2]*sizeof(double));
    }
    else
    {
        Iout=memcpy(Iout,mxGetData(prhs[5]),Osize[0]*Osize[1]*Osize[2]*sizeof(double));
    }
    
    for(q=0; q<Tsize[0]; q++)
    {
        t0=(int)TRI[q]-1;
        t1=(int)TRI[q+Tsize[0]]-1;
        t2=(int)TRI[q+Tsize[0]*2]-1;
        
        /* Vertices */
        p0[0]=UV[t0]-1; p0[1]=UV[t0+Vsize[0]]-1;
        p1[0]=UV[t1]-1; p1[1]=UV[t1+Vsize[0]]-1;
        p2[0]=UV[t2]-1; p2[1]=UV[t2+Vsize[0]]-1;
        
        /* Vertices2*/
        x0[0]=XY[t0]-1; x0[1]=XY[t0+Vsize[0]]-1;
        x1[0]=XY[t1]-1; x1[1]=XY[t1+Vsize[0]]-1;
        x2[0]=XY[t2]-1; x2[1]=XY[t2+Vsize[0]]-1;
        
        /*  Get bounding box (ROI) */
        boundmin[0]=(int)floor(min(min(p0[0],p1[0]),p2[0]));
        boundmin[1]=(int)floor(min(min(p0[1],p1[1]),p2[1]));
        
        boundmax[0]=(int) ceil(max(max(p0[0],p1[0]),p2[0]));
        boundmax[1]=(int) ceil(max(max(p0[1],p1[1]),p2[1]));
        
        boundmin[0]=max(boundmin[0],0);
        boundmin[1]=max(boundmin[1],0);
        
        boundmax[0]=min(boundmax[0],Osize[0]-1);
        boundmax[1]=min(boundmax[1],Osize[1]-1);
        
        /* Normalization factors */
        f12 = ( p1[1] - p2[1] ) * p0[0]  + (p2[0] - p1[0] ) * p0[1] + p1[0] * p2[1] - p2[0] *p1[1];
        f20 = ( p2[1] - p0[1] ) * p1[0]  + (p0[0] - p2[0] ) * p1[1] + p2[0] * p0[1] - p0[0] *p2[1];
        f01 = ( p0[1] - p1[1] ) * p2[0]  + (p1[0] - p0[0] ) * p2[1] + p0[0] * p1[1] - p1[0] *p0[1];
        
        /* Lambda Gradient */
        g12[0]=( p1[1] - p2[1] )/f12; g12[1] = (p2[0] - p1[0] )/f12;
        g20[0]=( p2[1] - p0[1] )/f20; g20[1] = (p0[0] - p2[0] )/f20;
        g01[0]=( p0[1] - p1[1] )/f01; g01[1] = (p1[0] - p0[0] )/f01;
        
        /* Center compensation */
        c12 = (p1[0] * p2[1] - p2[0] *p1[1])/f12;
        c20 = (p2[0] * p0[1] - p0[0] *p2[1])/f20;
        c01 = (p0[0] * p1[1] - p1[0] *p0[1])/f01;
        
        Lambdat[0]=g12[1]*boundmin[1]+c12+g12[0]*boundmin[0];;
        Lambdat[1]=g20[1]*boundmin[1]+c20+g20[0]*boundmin[0];;
        Lambdat[2]=g01[1]*boundmin[1]+c01+g01[0]*boundmin[0];;
        for(j=boundmin[1]; j<=boundmax[1]; j++)
        {
            Lambda[0] = Lambdat[0];
            Lambda[1] = Lambdat[1];
            Lambda[2] = Lambdat[2];
            for(i=boundmin[0]; i<=boundmax[0]; i++)
            {
                /* Check if voxel is inside the triangle */
                if((Lambda[0]>mval)&&(Lambda[0]<bval)&&(Lambda[1]>mval)&&(Lambda[1]<bval)&&(Lambda[2]>mval)&&(Lambda[2]<bval))
                {
                    posuv[0]=Lambda[0]*x0[0]+Lambda[1]*x1[0]+Lambda[2]*x2[0];
                    posuv[1]=Lambda[0]*x0[1]+Lambda[1]*x1[1]+Lambda[2]*x2[1];
                    
                    if(Osize[2]>1)
                    {
                        interpolate_2d_double_color(Ipixel,posuv[0], posuv[1], Isize, Iin,true,false);
                    }
                    else
                    {
                        Ipixel[0]=interpolate_2d_double_gray(posuv[0], posuv[1], Isize, Iin,true,false);
                    }
                    
                    for(rgb=0;rgb<Osize[2];rgb++)
                    {
                        Iout[i+j*Osize[0]+rgb*Oslice]=Ipixel[rgb];
                    }
                }
                Lambda[0] += g12[0];
                Lambda[1] += g20[0];
                Lambda[2] += g01[0];
            }
            Lambdat[0] += g12[1];
            Lambdat[1] += g20[1];
            Lambdat[2] += g01[1];
        }
    }
}