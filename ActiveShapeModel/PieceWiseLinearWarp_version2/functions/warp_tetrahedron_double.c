#include "mex.h"
#include "math.h"
#include "image_interpolation.h"
#include "multiple_os_thread.h"
#include "string.h"

/* The matlab mex function */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[] )
{
    double *Iout, *Iin, *XYZ, *UVW, *TETRA, *SizeO;
    float *Iinf, *Ioutf;
    const mwSize *Idims;
    const mwSize *Tdims;
    const mwSize *Vdims;
    
    int Isize[3] = {1,1,1};
    int Osize[3] = {1,1,1};
    int Tsize[2] = {1,1};
    int Vsize[2] = {1,1};
    
    /* loop variable */
    int q, i, j, k;
    
    /* Current Vertices indices*/
    int t0, t1, t2, t3;
    
    /* Index of current voxel/pixel */
    int index;
    
    /* Bounding box polygon */
    int boundmin[3];
    int boundmax[3];
    
    /* Vertices */
    double p0[3],p1[3], p2[3], p3[3];
    double x0[3],x1[3], x2[3], x3[3];
    
    /* Barycentric variables */
    double posuvw[3];
    
    /* Interpolation percentages */
    double Lambda[4],Lambdar[3],Lambdat[3];
    
    /* Transformation matrix */
    double T[9], Tinv[9], detT, detTi;
    
    /* Center compensation */
    double c1,c2,c3;
    
	/* Image class */
	int classdouble;
    
    /* Speed up */
    int checkair, ko;
	
	/* Check for proper number of arguments. */
	if(nrhs<5) {
		mexErrMsgTxt("Five inputs are required.");
	}
 
     /* function J=warp_triangle_double(I,xy,uv,tri,ImageSize) */
    if(mxIsDouble(prhs[0]))
	{
	    classdouble=true;
	}
	else if(mxIsSingle(prhs[0]))
	{ 
		classdouble=false;
	}
	else
	{
	    mexErrMsgTxt("Image must be single or double");
	}
	
	if(classdouble)
	{
		Iin=mxGetPr(prhs[0]);
	}
	else
	{
		Iinf=(float *)mxGetData(prhs[0]);
	}
	
    XYZ=mxGetPr(prhs[1]);
    UVW=mxGetPr(prhs[2]);
    TETRA=mxGetPr(prhs[3]);
    SizeO=mxGetPr(prhs[4]);
    
    /* Input Image size */
    Idims = mxGetDimensions(prhs[0]);
    Isize[0] = Idims[0];
    Isize[1] = Idims[1];
    Isize[2] = Idims[2];
    
    
    /* Input Number of Polygons */
    Tdims = mxGetDimensions(prhs[3]);
    Tsize[0] = Tdims[0];
    Tsize[1] = Tdims[1];
    
    /* Input Number of Polygons */
    Vdims = mxGetDimensions(prhs[1]);
    Vsize[0] = Vdims[0];
    Vsize[1] = Vdims[1];
    
    /* Output Image size */
    Osize[0]=(int)SizeO[0];
    Osize[1]=(int)SizeO[1];
    Osize[2]=(int)SizeO[2];
    
   
    /* Create empty array for output */
	if(classdouble)
	{
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
        
    }
	else
	{
		plhs[0] = mxCreateNumericArray(3, Osize, mxSINGLE_CLASS, mxREAL);
		Ioutf=(float *)mxGetData(plhs[0]);
        if(nrhs<6)
        {
            Ioutf=memset(Ioutf,0,Osize[0]*Osize[1]*Osize[2]*sizeof(float));
        }
        else
        {
            Ioutf=memcpy(Ioutf,mxGetData(prhs[5]),Osize[0]*Osize[1]*Osize[2]*sizeof(float));
        }
	}
    
    for(q=0; q<Tsize[0]; q++)
    {
        t0=(int)TETRA[q]-1;
        t1=(int)TETRA[q+Tsize[0]]-1;
        t2=(int)TETRA[q+Tsize[0]*2]-1;
        t3=(int)TETRA[q+Tsize[0]*3]-1;
        
        /* Vertices */
        p0[0]=UVW[t0]-1; p0[1]=UVW[t0+Vsize[0]]-1; p0[2]=UVW[t0+Vsize[0]*2]-1;
        p1[0]=UVW[t1]-1; p1[1]=UVW[t1+Vsize[0]]-1; p1[2]=UVW[t1+Vsize[0]*2]-1;
        p2[0]=UVW[t2]-1; p2[1]=UVW[t2+Vsize[0]]-1; p2[2]=UVW[t2+Vsize[0]*2]-1;
        p3[0]=UVW[t3]-1; p3[1]=UVW[t3+Vsize[0]]-1; p3[2]=UVW[t3+Vsize[0]*2]-1;
        
        /* Vertices2*/
        x0[0]=XYZ[t0]-1; x0[1]=XYZ[t0+Vsize[0]]-1; x0[2]=XYZ[t0+Vsize[0]*2]-1;
        x1[0]=XYZ[t1]-1; x1[1]=XYZ[t1+Vsize[0]]-1; x1[2]=XYZ[t1+Vsize[0]*2]-1;
        x2[0]=XYZ[t2]-1; x2[1]=XYZ[t2+Vsize[0]]-1; x2[2]=XYZ[t2+Vsize[0]*2]-1;
        x3[0]=XYZ[t3]-1; x3[1]=XYZ[t3+Vsize[0]]-1; x3[2]=XYZ[t3+Vsize[0]*2]-1;
        
        /*  Get bounding box (ROI) */
        boundmin[0]=(int)floor(min(min(min(p0[0],p1[0]),p2[0]),p3[0]));
        boundmin[1]=(int)floor(min(min(min(p0[1],p1[1]),p2[1]),p3[1]));
        boundmin[2]=(int)floor(min(min(min(p0[2],p1[2]),p2[2]),p3[2]));
        
        boundmax[0]=(int) ceil(max(max(max(p0[0],p1[0]),p2[0]),p3[0]));
        boundmax[1]=(int) ceil(max(max(max(p0[1],p1[1]),p2[1]),p3[1]));
        boundmax[2]=(int) ceil(max(max(max(p0[2],p1[2]),p2[2]),p3[2]));
        
        boundmin[0]=max(boundmin[0],0);
        boundmin[1]=max(boundmin[1],0);
        boundmin[2]=max(boundmin[2],0);
        
        boundmax[0]=min(boundmax[0],Osize[0]-1);
        boundmax[1]=min(boundmax[1],Osize[1]-1);
        boundmax[2]=min(boundmax[2],Osize[2]-1);
        
        /*  Define matrix T */
        T[0]=p0[0]-p3[0]; T[3]=p1[0]-p3[0]; T[6]= p2[0]-p3[0];
        T[1]=p0[1]-p3[1]; T[4]=p1[1]-p3[1]; T[7]= p2[1]-p3[1];
        T[2]=p0[2]-p3[2]; T[5]=p1[2]-p3[2]; T[8]= p2[2]-p3[2];
        
        /* matrix inverse of T */
        detT=T[0]*(T[4]*T[8]-T[7]*T[5])+ T[3]*(T[7]*T[2]-T[8]*T[1])+ T[6]*(T[1]*T[5]-T[4]*T[2]);
        
        detTi=1/detT;
        Tinv[0]=detTi*(T[4]*T[8]-T[5]*T[7]); Tinv[3]=detTi*(T[5]*T[6]-T[3]*T[8]); Tinv[6]=detTi*(T[3]*T[7]-T[4]*T[6]);
        Tinv[1]=detTi*(T[7]*T[2]-T[8]*T[1]); Tinv[4]=detTi*(T[8]*T[0]-T[6]*T[2]); Tinv[7]=detTi*(T[6]*T[1]-T[7]*T[0]);
        Tinv[2]=detTi*(T[1]*T[5]-T[2]*T[4]); Tinv[5]=detTi*(T[2]*T[3]-T[0]*T[5]); Tinv[8]=detTi*(T[0]*T[4]-T[1]*T[3]);
        
        /* Center compensation */
        c1=-Tinv[0]*p3[0]-Tinv[3]*p3[1]-Tinv[6]*p3[2];
        c2=-Tinv[1]*p3[0]-Tinv[4]*p3[1]-Tinv[7]*p3[2];
        c3=-Tinv[2]*p3[0]-Tinv[5]*p3[1]-Tinv[8]*p3[2];
                        
        Lambdar[0]=Tinv[6]*boundmin[2]+c1+Tinv[0]*boundmin[0]+Tinv[3]*boundmin[1];
        Lambdar[1]=Tinv[7]*boundmin[2]+c2+Tinv[1]*boundmin[0]+Tinv[4]*boundmin[1];
        Lambdar[2]=Tinv[8]*boundmin[2]+c3+Tinv[2]*boundmin[0]+Tinv[5]*boundmin[1];
        
        for(k=boundmin[2]; k<=boundmax[2]; k++)
        {
            Lambdat[0]=Lambdar[0];
            Lambdat[1]=Lambdar[1];
            Lambdat[2]=Lambdar[2];
            ko = k*Osize[0]*Osize[1];
            for(j=boundmin[1]; j<=boundmax[1]; j++)
            {
                Lambda[0]=Lambdat[0];
                Lambda[1]=Lambdat[1];
                Lambda[2]=Lambdat[2];
                checkair=false;
                for(i=boundmin[0]; i<=boundmax[0]; i++)
                {
                    /* Check if voxel is inside the tetrahedron */
                    if((Lambda[0]>-0.0001)&&(Lambda[0]<1.0001)&&(Lambda[1]>-0.0001)&&(Lambda[1]<1.0001)&&(Lambda[2]>-0.0001)&&(Lambda[2]<1.0001))
                    {
                        Lambda[3]=1-Lambda[0]-Lambda[1]-Lambda[2];
                        if((Lambda[3]>-0.0001)&&(Lambda[3]<1.0001))
                        {
                            posuvw[0]=Lambda[0]*x0[0]+Lambda[1]*x1[0]+Lambda[2]*x2[0]+Lambda[3]*x3[0];
                            posuvw[1]=Lambda[0]*x0[1]+Lambda[1]*x1[1]+Lambda[2]*x2[1]+Lambda[3]*x3[1];
                            posuvw[2]=Lambda[0]*x0[2]+Lambda[1]*x1[2]+Lambda[2]*x2[2]+Lambda[3]*x3[2];

                            index=i+j*Osize[0]+ko;
                            if(classdouble)
                            {
                                Iout[index]=interpolate_3d_double_gray(posuvw[0], posuvw[1], posuvw[2], Isize, Iin,false,false);
                            }
                            else
                            {
                                Ioutf[index]=interpolate_3d_float_gray((float)posuvw[0], (float)posuvw[1], (float)posuvw[2], Isize, Iinf,false,false);
                            }
                            checkair=true;
                        }
                    }
                    else if(checkair)
                    {
                        break;
                    }
                    Lambda[0]+=Tinv[0];
                    Lambda[1]+=Tinv[1];
                    Lambda[2]+=Tinv[2];
                }
                Lambdat[0]+=Tinv[3];
                Lambdat[1]+=Tinv[4];
                Lambdat[2]+=Tinv[5];
            }
            Lambdar[0]+=Tinv[6];
            Lambdar[1]+=Tinv[7];
            Lambdar[2]+=Tinv[8];
        }
    }
}