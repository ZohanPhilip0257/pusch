#include<stdio.h>
#include<math.h>
#include"wnDataTypes.h"
#define PI                  ( (wnFlt) 3.14159265358979 )

wnFlt sinLUT[129] = {0.000000, 0.012272, 0.024541, 0.036807, 0.049068,\
                    0.061321, 0.073565, 0.085797, 0.098017, 0.110222,\
                    0.122411, 0.134581, 0.146730, 0.158858, 0.170962,\
                    0.183040, 0.195090, 0.207111, 0.219101, 0.231058,\
                    0.242980, 0.254866, 0.266713, 0.278520, 0.290285,\
                    0.302006, 0.313682, 0.325310, 0.336890, 0.348419,\
                    0.359895, 0.371317, 0.382683, 0.393992, 0.405241,\
                    0.416430, 0.427555, 0.438616, 0.449611, 0.460539,\
                    0.471397, 0.482184, 0.492898, 0.503538, 0.514103,\
                    0.524590, 0.534998, 0.545325, 0.555570, 0.565732,\
                    0.575808, 0.585798, 0.595699, 0.605511, 0.615232,\
                    0.624859, 0.634393, 0.643832, 0.653173, 0.662416,\
                    0.671559, 0.680601, 0.689541, 0.698376, 0.707107,\
                    0.715731, 0.724247, 0.732654, 0.740951, 0.749136,\
                    0.757209, 0.765167, 0.773010, 0.780737, 0.788346,\
                    0.795837, 0.803208, 0.810457, 0.817585, 0.824589,\
                    0.831470, 0.838225, 0.844854, 0.851355, 0.857729,\
                    0.863973, 0.870087, 0.876070, 0.881921, 0.887640,\
                    0.893224, 0.898674, 0.903989, 0.909168, 0.914210,\
                    0.919114, 0.923880, 0.928506, 0.932993, 0.937339,\
                    0.941544, 0.945607, 0.949528, 0.953306, 0.956940,\
                    0.960431, 0.963776, 0.966976, 0.970031, 0.972940,\
                    0.975702, 0.978317, 0.980785, 0.983105, 0.985278,\
                    0.987301, 0.989177, 0.990903, 0.992480, 0.993907,\
                    0.995185, 0.996313, 0.997290, 0.998118, 0.998795,\
                    0.999322, 0.999699, 0.999925, 1.000000};

wnFlt wnSinRad(double inp)
{
    wnFlt delX = PI*0.00390625;//pi/256
    wnUInt16 secNum;
    wnFlt x1, x2, y1, y2;
    wnFlt mulFac = 1;
    
    if(inp>2*PI || inp<-2*PI)
    {
        inp = fmod(inp, 2*PI);
    }

    /*while(inp>2*PI)
    {
        inp = inp-2*PI;
    }
    while(inp<-2*PI)
    {
        inp = inp+2*PI;
    }//*/

    if(inp<0)
    {
        inp = -1*inp;
        mulFac = -1;
    }

    if(inp>=0 && inp<=PI*0.5)
    {
        if(inp==PI*0.5)
        {
            return sinLUT[128]*mulFac;
        }
        else
        {
            secNum = (inp/delX);
            x1 = secNum*delX;
            x2 = (secNum+1)*delX;
            y1 = sinLUT[secNum];
            y2 = sinLUT[secNum+1];
            return (y1+(y2-y1)/(x2-x1)*(inp-x1))*mulFac;
        }
    }
    else if(inp>PI*0.5 && inp<=PI)
    {
        inp = PI-inp;
        if(inp==PI*0.5)
        {
            return sinLUT[128]*mulFac;
        }
        else
        {
            secNum = (inp/delX);
            x1 = secNum*delX;
            x2 = (secNum+1)*delX;
            y1 = sinLUT[secNum];
            y2 = sinLUT[secNum+1];
            return (y1+(y2-y1)/(x2-x1)*(inp-x1))*mulFac;

        }
    }
    else if(inp>PI && inp<=PI*1.5)
    {
        inp = inp-PI;
        if(inp==PI*0.5)
        {
            return sinLUT[128]*mulFac*-1;
        }
        else
        {
            secNum = (inp/delX);
            x1 = secNum*delX;
            x2 = (secNum+1)*delX;
            y1 = sinLUT[secNum];
            y2 = sinLUT[secNum+1];
            return (-1*(y1+(y2-y1)/(x2-x1)*(inp-x1)))*mulFac;
        }
    }
    else
    {
        if(inp==PI*0.5)
        {
            return sinLUT[128]*mulFac*-1;
        }
        else
        {
            inp = 2*PI - inp;
            secNum = (inp/delX);
            x1 = secNum*delX;
            x2 = (secNum+1)*delX;
            y1 = sinLUT[secNum];
            y2 = sinLUT[secNum+1];
            return (-1*(y1+(y2-y1)/(x2-x1)*(inp-x1)))*mulFac;
        }
    }
}

wnFlt wnCosRad(double inp)
{
    return wnSinRad(inp+PI*0.5);
}