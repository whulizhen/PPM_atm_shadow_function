//============================================================================
//
//  This file is part of GFC, the GNSS FOUNDATION CLASS.
//
//  The GFC is free software; you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published
//  by the Free Software Foundation; either version 3.0 of the License, or
//  any later version.
//
//  The GFC is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with GFC; if not, write to the Free Software Foundation,
//  Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110, USA
//
//  Copyright 2017, lizhen
//
//============================================================================


#include "GVector.hpp"
#include <cmath>
#include <complex>
using namespace gfc;
using namespace std;


    // solve cubic equation x^3 + a*x^2 + b*x + c
    // x - array of size 3
    // In case 3 real roots: => x[0], x[1], x[2], return 3
    //         2 real roots: x[0], x[1],          return 2
    //         1 real root : x[0], x[1] ± i*x[2], return 1
    unsigned int solveP3(double *x,double a,double b,double c)
    {
        double a2 = a*a;
        double q  = (a2 - 3*b)/9;
        double r  = (a*(2*a2-9*b) + 27*c)/54;
        double r2 = r*r;
        double q3 = q*q*q;
        double A,B;
        double M_2PI = M_PI*2;
        double eps = 1.0e-12;
        if(r2<q3)
        {
            double t=r/sqrt(q3);
            if( t<-1) t=-1;
            if( t> 1) t= 1;
            t=acos(t);
            a/=3; q=-2*sqrt(q);
            x[0]=q*cos(t/3)-a;
            x[1]=q*cos((t+M_2PI)/3)-a;
            x[2]=q*cos((t-M_2PI)/3)-a;
            return 3;
        }
        else
        {
            A =-pow(fabs(r)+sqrt(r2-q3),1./3);
            if( r<0 ) A=-A;
            B = (0==A ? 0 : q/A);

            a/=3;
            x[0] =(A+B)-a;
            x[1] =-0.5*(A+B)-a;
            x[2] = 0.5*sqrt(3.)*(A-B);
            if(fabs(x[2])<eps) { x[2]=x[1]; return 2; }

            return 1;
        }
    }

// solve quartic equation x^4 + a*x^3 + b*x^2 + c*x + d
    unsigned int solve_quartic(double a, double b, double c, double d,double x[4] )
    {
        double eps = 1e-5;
        int num_real = 0;
        std::complex<double> retval[4];

        // https://math.stackexchange.com/questions/785/is-there-a-general-formula-for-solving-4th-degree-equations-quartic
        double A = b - 3.0/8.0*a*a;
        double B = c - a*b/2.0 + a*a*a/8.0;
        double C = d - a*c/4.0 + a*a*b/16.0 - 3*a*a*a*a/256;

        double coef3[3]={0.0};
        coef3[0] = -A/2.0;
        coef3[1] = -C;
        coef3[2] = (4*A*C-B*B)/8;
        double x3[3];
        unsigned int iZeroes = solveP3(x3, coef3[0], coef3[1], coef3[2]);
        double s = x3[0];
        // The essence - choosing Y with maximal absolute value.
        if(iZeroes != 1)
        {
            if( fabs(x3[1]) > fabs(s)) s = x3[1];
            if( fabs(x3[2]) > fabs(s)) s = x3[2];
        }


        // what is the solution ? A= 2s, i.e. B==0
        // y^4 + Ay^2 + C = 0
        if( fabs(B)< eps)
        {
            double D = A*A -4*C;
            complex<double> x1 = (-A + sqrt(D))*0.5;
            complex<double> x2 = (-A - sqrt(D))*0.5;

            retval[0] = sqrt(x1)- a/4;
            retval[1] = -sqrt(x1)- a/4;
            retval[2] =  sqrt(x2) - a/4;
            retval[3] = -sqrt(x2) -a/4;
        }
        else
        {
            complex<double> sqD1 = sqrt(2*s -A);
            complex<double> sqD2 = sqrt(-2*s - A + 2*B/sqD1);
            complex<double> sqD3 = sqrt(-2*s - A - 2*B/sqD1);

            retval[0] = -0.5*sqD1 + 0.5*sqD2 - a/4.0;
            retval[1] = -0.5*sqD1 - 0.5*sqD2 - a/4.0;
            retval[2] =  0.5*sqD1 + 0.5*sqD3 - a/4.0;
            retval[3] =  0.5*sqD1 - 0.5*sqD3 - a/4.0;
        }

        for(int i = 0 ; i< 4; i++)
        {
            if( fabs(retval[i].imag())< eps )
            {
                x[num_real] = retval[i].real();
                num_real++;
            }
        }

       return num_real;

   }



      double area_hyperbola(double Q1[2], double Q2[2], double Os[2],double Oe[2], double Rs,double a, double b, bool in_out, bool x_axis)
      {
        double area = 0.0;
        double as[2] = {Q1[0]-Os[0],Q1[1] -Os[1] };
        double bs[2] = {Q2[0]-Os[0],Q2[1] -Os[1] };
        
        double ae[2] = {Q1[0]-Oe[0],Q1[1] -Oe[1]  };
        double be[2] = {Q2[0]-Oe[0],Q2[1] -Oe[1]  };
        
        
        double TQ1Q2Os = 0.5*fabs(as[0]*bs[1] - bs[0]*as[1]);
        double TQ1Q2Oe = 0.5*fabs(ae[0]*be[1] - be[0]*ae[1]);
        
        double cts = (as[0]*bs[0] + as[1]*bs[1] )/sqrt( (as[0]*as[0]+as[1]*as[1]) * (bs[0]*bs[0]+bs[1]*bs[1]) );
        double SQ1Q2Os = 0.5 * acos(cts)*Rs*Rs;
        
        //calculate the area of hyperbolic secttion
        
        if(x_axis == false) //// 主轴是y, 整个双曲线旋转90度，使得x为主轴
        {
            // rotation transformation
            // x' = y
            // y' =-x
            double tmp = 0.0;
            tmp = Q1[0]; Q1[0] = Q1[1]; Q1[1] = -tmp;
            tmp = Q2[0]; Q2[0] = Q2[1]; Q2[1] = -tmp;
            tmp = a;
            a = b;
            b = tmp;
        }
        
        //ref: http://www.robertobigoni.eu/Matematica/Conics/segmentHyp/segmentHyp.html
        double s1 = 0.0, s2 =0.0;
        
        double xQ1 = fabs(Q1[0]);
        double xQ2 = fabs(Q2[0]);
        
        s1 = b*( xQ1*sqrt(  (xQ1*xQ1/a/a) -1.0 ) - a*acosh(xQ1/a) );
        
        s2 = b*( xQ2*sqrt(  (xQ2*xQ2/a/a) -1.0 ) - a*acosh(xQ2/a) );
        
        if(isnan(s1) || isnan(s2) )
        {
            int testc = 0;
        }
        
        double ss = s1<=s2?s1:s2;
        double sl = s1>=s2?s1:s2;
        double S2 =0.0;
        if( Q1[1]*Q2[1] < 0.0 )  //2个交点在y轴不同侧
        {
            double x_m = 0.0;
            
            if(fabs(Q2[0] - Q1[0])>1.0E-10)
            {
                double k = ( Q2[1] - Q1[1] )/( Q2[0] - Q1[0]);
                double m = Q1[1] - k *Q1[0];
                x_m = - m/k;
            }
            else
            {
                x_m = Q1[0];
            }
            
            if(Q1[0]>Q2[0])
            {
                S2 = ss + (sl - ss)/2.0 - 0.5*fabs((x_m - Q1[0])*Q1[1]) + 0.5*fabs((Q2[0]- x_m)*Q2[1]) ;
            }
            else
            {
                S2 = ss + (sl - ss)/2.0 - 0.5*fabs((Q2[0]- x_m)*Q2[1]) + 0.5*fabs((x_m - Q1[0])*Q1[1]);
            }
            
            
        }
        else  //2 个交点在同一侧
        {
            double s_trapizium = 0.5*fabs( Q1[0] - Q2[0] )* (fabs(Q1[1]) + fabs(Q2[1]) );
            S2 = (sl - ss - 2*s_trapizium )/2.0;
        }
        
        
        if(in_out == true) // the centre of the sun is inside the ellipse
        {
            area = SQ1Q2Os - TQ1Q2Os - S2;
        }
        else // // the centre of the sun is outside the ellipse
        {
            double area_shadow =  SQ1Q2Os - TQ1Q2Os + S2;
            area = 3.14159265357*Rs*Rs -area_shadow;
        }
        
        return area;

      }

      double area_ellispe(double Q1[2], double Q2[2], double Os[2], double Rs,double a, double b, bool in_out)
      {
          double area = 0.0;

        double as[2] = {Q1[0]-Os[0],Q1[1] -Os[1] };
        double bs[2] = {Q2[0]-Os[0],Q2[1] -Os[1] };

        double ae[2] = {Q1[0],Q1[1]  };
        double be[2] = {Q2[0],Q2[1]  };


        double TQ1Q2Os = 0.5*fabs(as[0]*bs[1] - bs[0]*as[1]);
        double TQ1Q2Oe = 0.5*fabs(ae[0]*be[1] - be[0]*ae[1]);

        double cts = (as[0]*bs[0] + as[1]*bs[1] )/sqrt( (as[0]*as[0]+as[1]*as[1]) * (bs[0]*bs[0]+bs[1]*bs[1]) );
        double SQ1Q2Os = 0.5 * acos(cts)*Rs*Rs;

        //
        //double cte = (ae[0]*be[0] + ae[1]*be[1] )/sqrt( (ae[0]*ae[0]+ae[1]*ae[1]) * (be[0]*be[0]+be[1]*be[1]) );
        // elliptical sector
        //double SQ1Q2Oe = 0.5*a*b*acos(cte);

        double SQ1Q2Oe = 0.0;
        double aa =  a> b? a: b;
        double q1 =0.0, q2 =0.0;
        if(ae[0]>0 && ae[1]>0)  // the first quadrant
        {
            q1 = atan(ae[1]/ ae[0]);
        }
        else if(ae[0]<0 && ae[1]>0)
        {
            q1 = 3.14159265357 - atan(-ae[1]/ae[0]);
        }
        else if(ae[0]<0 && ae[1]<0)
        {
            q1 = 3.14159265357 + atan(ae[1]/ ae[0]);
        }
        else if( ae[0]>0 && ae[1]<0 )
        {
            q1 = 3.14159265357*2 - atan(-ae[1]/ae[0]);
        }

        if(be[0]>0 && be[1]>0)  // the first quadrant
        {
            q2 = atan(be[1]/be[0]);
        }
        else if(be[0]<0 && be[1]>0)
        {
            q2 = 3.14159265357 - atan(-be[1]/be[0]);
        }
        else if(be[0]<0 && be[1]<0)
        {
            q2 = 3.14159265357 + atan(be[1]/be[0]);
        }
        else if( be[0]>0 && be[1]<0 )
        {
            q2 = 3.14159265357*2 - atan(-be[1]/be[0]);
        }

        double s_a = a*b/2.0*( q1 - atan2( (b-a)*sin(2.0*q1), (a+b) + (b-a)*cos(2*q1)) );
        double s_b = a*b/2.0*( q2 - atan2( (b-a)*sin(2.0*q2), (a+b) + (b-a)*cos(2*q2)) );

        SQ1Q2Oe = fabs(s_b - s_a);

        double S1 = SQ1Q2Oe - TQ1Q2Oe;

        if(in_out == true) // the centre of the sun is inside the ellipse
        {
            area = SQ1Q2Os - TQ1Q2Os - S1;
        }
        else // // the centre of the sun is outside the ellipse
        {
            double area_shadow =  SQ1Q2Os - TQ1Q2Os + S1;
            area = 3.14159265357*Rs*Rs -area_shadow;
        }

        return area;
      }


int myperspectiveProjection(double a, double b, GVector& sunpos_ecef, GVector& satpos_ecef, double& r_solar, double& area_bright, double& dis_boundary, double& dis_circle)
   {
       int state = -1;
        double Rs = 695700.0; //km
        double ab = a*b;
        double a2 = a*a;
        double b2 = b*b;
        double ab2 = ab*ab;

        GVector r = satpos_ecef;
        GVector rs = sunpos_ecef;

        double t =0.0, t1 = 0.0, t2 =0.0, dis =0, s1 = 0, s2=0,ds1 = 0.0, ds2 =0.0;
        //A: first test if the satellite is in the front of earth,
        // if it is in the front of earth, it is always full phase
        // otherwise,using the photogrammetry method
        // A= diag{1/a2,1/a2, 1/b2 }, A^{-1} = diag{a2, a2, b2}

        GVector o, u, v, n;
        double f = 1000.0;
        n = r - rs;
        double dis_sat_earth = r.norm();
        double dis_sat_sun = n.norm();
        n.normalise();

        double nAin = n.x*n.x*a2 + n.y*n.y*a2 + + n.z*n.z*b2;
        double rtn = dotproduct(r, n);
        double ntn = dotproduct(n, n);


        s1 = sqrt(1.0/nAin);
        s2 = - s1;

        t1 = (rtn - 1.0/s1)/ntn;
        t2 = (rtn - 1.0/s2)/ntn;

        GVector xs1 = r - t1*n;
        GVector xs2 = r - t2*n;

        ds1 = (xs1-rs).norm();
        ds2 = (xs2-rs).norm();


        t = (ds1 <= ds2) ? ds1:ds2;
        ds2 = (ds1 >= ds2)? ds1:ds2;
        ds1 = t;

        double nAn = n.x*n.x/a2 + n.y*n.y/a2 + n.z*n.z/b2;
        double rsAn = n.x*rs.x/a2 + n.y*rs.y/a2 + n.z*rs.z/b2;
        double rsArs = rs.x*rs.x/a2 + rs.y*rs.y/a2 + rs.z*rs.z/b2;
        double Delta = rsAn*rsAn - nAn*(rsArs-1.0);

        if(Delta > 0) // sun-sat line intersects the Earth ellipsoid
        {
            t1 = (-2.0*rsAn + sqrt(Delta))/2.0/nAn;
            t2 = (-2.0*rsAn - sqrt(Delta))/2.0/nAn;

            ds1 = (t1 <= t2) ? t1:t2;
            //ds2 = (t1 <= t2) ? t2:t1;
        }
        else
        {
            // normal vector to the plane sat,sun and Earth
            GVector p = crossproduct(r, n);
            // normal at the ellipsoid that is perpendicular to n
            GVector ns = crossproduct(n, p);
            ns.normalise();

            t = ns.x*a2*ns.x + ns.y*a2*ns.y + ns.z*b2*ns.z ;
            double lam1 = sqrt( 1.0/ t );
            double lam2 = -lam1;

            s1 = dotproduct(ns, rs) - lam1*t;
            s2 = dotproduct(ns, rs) - lam2*t;

            double lam = fabs(s1)<fabs(s2)? lam1: lam2;
            dis = lam*(n.x*ns.x*a2 + n.y*ns.y*a2 + n.z*ns.z*b2) - dotproduct(n, rs);

            ds1 = dis;
        }

        if(dis_sat_sun < ds1)  // full phase
        {
            state = 0;
            //导致后面的 area_bright, dis_boundary, and dis_circle 都没有计算

            return state;
        }
        else if(dis_sat_sun >= ds2)  // ellipse
        {
            state = 2;
        }
        else if( dis_sat_sun >= ds1 && dis_sat_sun <= ds2 ) // hyperbola
        {
             state = 1;  // hyperbola
        }

        //B: build the ISF coordinate system
        o = r - f*n; // the origin of the photo coordinate system (PSC)

        u = n - 1.0/rtn*r;

        u.normalise();
        v = crossproduct(n, u);

        // the projection of earth's centre PEC in ECEF
        double t_dis = sqrt( (f*dis_sat_earth/rtn)*(f*dis_sat_earth/rtn) - f*f );
        GVector PEC = t_dis*u; // vector from PSC to PEC
        GVector PEC_isf;
        // convert PEC into ISF
        PEC_isf.x = dotproduct(u, PEC);
        PEC_isf.y = dotproduct(v, PEC);
        PEC_isf.z = dotproduct(n, PEC);

        t = (r.x*r.x/a2 + r.y*r.y/a2 + r.z*r.z/b2 - 1.0);

        double M[3][3]={0.0};
        M[0][0] = r.x*r.x/a2/a2 - t/a2;
        M[0][1] = r.x*r.y/a2/a2;
        M[0][2] = r.x*r.z/a2/b2;
        M[1][0] = M[0][1];
        M[1][1] = r.y*r.y/a2/a2 - t/a2;
        M[1][2] = r.y*r.z/a2/b2;
        M[2][0] = M[0][2];
        M[2][1] = M[1][2];
        M[2][2] = r.z*r.z/b2/b2 - t/b2;


        //D: calculate K
        double K[6] = {0.0};

        K[0] = (u.x*M[0][0] + u.y*M[1][0] + u.z*M[2][0])*u.x
             + (u.x*M[0][1] + u.y*M[1][1] + u.z*M[2][1])*u.y
             + (u.x*M[0][2] + u.y*M[1][2] + u.z*M[2][2])*u.z;

        K[1] = (u.x*M[0][0] + u.y*M[1][0] + u.z*M[2][0])*v.x
            + (u.x*M[0][1] + u.y*M[1][1] + u.z*M[2][1])*v.y
            + (u.x*M[0][2] + u.y*M[1][2] + u.z*M[2][2])*v.z;

        K[2] = (v.x*M[0][0] + v.y*M[1][0] + v.z*M[2][0])*v.x
            + (v.x*M[0][1] + v.y*M[1][1] + v.z*M[2][1])*v.y
            + (v.x*M[0][2] + v.y*M[1][2] + v.z*M[2][2])*v.z;

        K[3] = -((u.x*M[0][0] + u.y*M[1][0] + u.z*M[2][0])*n.x
            + (u.x*M[0][1] + u.y*M[1][1] + u.z*M[2][1])*n.y
            + (u.x*M[0][2] + u.y*M[1][2] + u.z*M[2][2])*n.z)*f*2.0;

        K[4] = -((v.x*M[0][0] + v.y*M[1][0] + v.z*M[2][0])*n.x
            + (v.x*M[0][1] + v.y*M[1][1] + v.z*M[2][1])*n.y
            + (v.x*M[0][2] + v.y*M[1][2] + v.z*M[2][2])*n.z)*f*2.0;

        K[5] = ((n.x*M[0][0] + n.y*M[1][0] + n.z*M[2][0])*n.x
            + (n.x*M[0][1] + n.y*M[1][1] + n.z*M[2][1])*n.y
            + (n.x*M[0][2] + n.y*M[1][2] + n.z*M[2][2])*n.z)*f*f;


        // just make the numbers larger for visualisation
        for( int i =0; i< 6; i++ )
        {
            K[i] = K[i]*1.0E6;
        }

        //E: calculate the eigen value and eigen vector of matrix B
        t =  sqrt( (K[0]+K[2])*(K[0]+K[2]) - 4.0*(K[0]*K[2] - K[1]*K[1]) ) ;
        double lambda1 = (K[0]+K[2]+t)/2.0;
        double lambda2 = (K[0]+K[2]-t)/2.0;

        double  r1[2]={0.0,1.0},r2[2]={1.0,0.0};
        // http://www.math.harvard.edu/archive/21b_fall_04/exhibits/2dmatrices/index.html
        if( fabs(K[1]) < 1.0E-12)
        {
            r1[0] = 1;
            r1[1] = 0;

            r2[0] = 0;
            r2[1] = 1;
        }
        else
        {
            r1[0] = lambda1 - K[2];
            r1[1] = K[1];
            r2[0] = lambda2 - K[2];
            r2[1] = K[1];
        }

        //get the unit vector
        t = sqrt(r1[0]*r1[0]+r1[1]*r1[1]);
        r1[0] = r1[0]/t;r1[1] = r1[1]/t;

        t= sqrt(r2[0]*r2[0]+r2[1]*r2[1]);
        r2[0] = r2[0]/t;r2[1] = r2[1]/t;

        double OM = (K[2]*K[3]*K[3] - 2.0*K[1]*K[3]*K[4] + K[0]*K[4]*K[4] )/(4.0*(K[0]*K[2] - K[1]*K[1]));
        // the translation parameters
        double tx = (r1[0]*K[3] + r1[1]*K[4])/lambda1;
        double ty = (r2[0]*K[3] + r2[1]*K[4])/lambda2;

        // semi-major axis and the semi-minor axis
        double AA =  ( OM - K[5] )/lambda1;
        double BB =  ( OM - K[5] )/lambda2;
        double R0 = f*Rs/dis_sat_sun;
        r_solar = R0;

        //F: Solve the quartic in the translated and rotated frame
        double A = lambda1*(R0 - 0.5*tx )*(R0 - 0.5*tx) + 0.25*lambda2*ty*ty - OM + K[5];
        double B = 2.0*lambda2*R0*ty;
        double C = lambda1*(0.5*tx*tx - 2*R0*R0) + lambda2*(0.5*ty*ty + 4*R0*R0) - 2*OM + 2*K[5];
        double D = 2*lambda2*R0*ty;
        double E = lambda1*(R0+0.5*tx)*(R0+0.5*tx) + 0.25*lambda2*ty*ty - OM + K[5];

        double ce[4]={B/A,C/A,D/A,E/A};

        double XX[4]={0.0}; // the solution of quartic equation
        int num_of_solution = 0;


        num_of_solution = solve_quartic(ce[0], ce[1], ce[2], ce[3], XX);

        if( num_of_solution == 3 || num_of_solution == 4)
        {
            printf("WARNING: shadowfactor, %d intersections!\n", num_of_solution);
            exit(0);
        }

        //calculate the sun's projection centre Os (PSC) in the transformed frame
        double Os[2] = {0.5*tx , 0.5*ty };
        //double Os[2] = {0.0,0.0};
        double PEC_new[2] = {0.0,0.0};
        // the projection of the earth's centre Pe in the transformed frame is obtained by converting PEC_isf into the rotated and translated frame
        PEC_new[0] = PEC_isf.x*r1[0] + r1[1]*PEC_isf.y + 0.5*tx;
        PEC_new[1] = PEC_isf.x*r2[0] + r2[1]*PEC_isf.y + 0.5*ty;

        //printf("PEC_new: %f %f\n",PEC_new[0], PEC_new[1]);

        // calculate the coordinates of the intersections betweent the circle and the conical curve in transformed frame
        // in the new frame, the origin is at the centre of the projection of the Earth
        double Q1[2]={ (1.0 - XX[0]*XX[0])/(1.0 + XX[0]*XX[0])*R0 + 0.5*tx , 2*XX[0]/(1.0 + XX[0]*XX[0])*R0 + 0.5*ty};
        double Q2[2]={ (1.0 - XX[1]*XX[1])/(1.0 + XX[1]*XX[1])*R0 + 0.5*tx , 2*XX[1]/(1.0 + XX[1]*XX[1])*R0 + 0.5*ty};

        //figure out the intersection between the line from PEC_new to Os and the conical curve.
        // parameter equation of line from Oe to Os in transformed frame, x = td
        double PEC_Os_len = sqrt( (Os[0]-PEC_new[0])*(Os[0]-PEC_new[0]) + (Os[1]-PEC_new[1])*(Os[1]-PEC_new[1]));
        double d[2] = { (Os[0]-PEC_new[0])/PEC_Os_len,(Os[1]-PEC_new[1])/PEC_Os_len };
        double aa = lambda1*d[0]*d[0] + lambda2*d[1]*d[1];
        double bb = 2*(lambda1*PEC_new[0]*d[0] + lambda2*PEC_new[1]*d[1]);
        double cc = lambda1*PEC_new[0]*PEC_new[0] + lambda2*PEC_new[1]*PEC_new[1] + K[5] - OM;
        t1 = (-bb + sqrt(bb*bb - 4.0*aa*cc ))/aa/2.0;
        t2 = (-bb - sqrt(bb*bb - 4.0*aa*cc ))/aa/2.0;

        if(t1*t2<0.0)
        {
            t = t1>0?t1:t2;
        }
        else
        {
           t =  t1<=t2?t1:t2;
        }

        dis_boundary = t;

        //boundary_intersection[0] = PEC_new[0] + t*d[0];
        //boundary_intersection[1] = PEC_new[1] + t*d[1];

        // //figure out the intersection between the line from PEC_new to Os and the solar circle
        aa = d[0]*d[0] + d[1]*d[1];
        bb = 2*(PEC_new[0]*d[0]+PEC_new[1]*d[1]) - (d[0]*tx + d[1]*ty );
        cc = PEC_new[0]*PEC_new[0]+PEC_new[1]*PEC_new[1] - (PEC_new[0]*tx + PEC_new[1]*ty ) + 0.25*(tx*tx+ty*ty) - R0*R0;

        t1 = (-bb + sqrt(bb*bb - 4.0*aa*cc ))/aa/2.0;
        t2 = (-bb - sqrt(bb*bb - 4.0*aa*cc ))/aa/2.0;
        // t should be the smaller one, which means closer to the solid Earth
        if(t1*t2<0.0)
        {
            t = t1>0?t1:t2;
        }
        else
        {
            t =  t1<=t2?t1:t2;
        }

        dis_circle = t;


        t = K[0]*K[2] - K[1]*K[1];
        if( t > 0.0) // ellipse
        {
            bool in_out = false;

            if( OM/(OM-K[5]) <= 1.0 )
            {
                in_out = true;
            }
            else
            {
                in_out = false;
            }

            if( in_out == true && num_of_solution < 2 )  //umbra
            {
                state = -1;
                return state;
            }
            if( in_out == false && num_of_solution < 2 )  //full phase
            {
                state = 0;
                return state;
            }

            // penumbra
            area_bright = area_ellispe(Q1, Q2, Os, R0, sqrt(AA), sqrt(BB), in_out);

        }
        else if (t < 0.0) // hyperbola
        {
            bool in_out = false;

            if( OM/(OM-K[5]) <= 1.0 )
            {
                in_out = false;
            }
            else
            {
                in_out = true;
            }

            if( in_out == true && num_of_solution <2)
            {
                state = -1; // umbra
                return state;
            }

            if( in_out == false && num_of_solution <2)
            {
                state = 0; // full phase
                return state;
            }

            // penumbra
            double a =0, b =0;
            bool x_axis = true;
            if(AA > 0)
            {
                a =sqrt(AA);
            }
            else
            {
                a = sqrt(-AA);
                x_axis = false;
            }
            if(BB > 0)
            {
                b =sqrt(BB);
                x_axis =false;
            }
            else
            {
                b = sqrt(-BB);
            }

            area_bright = area_hyperbola( Q1, Q2, Os,PEC_new, R0, a, b, in_out, x_axis);

        }

        return state;

   }




// the main entry of the shadow function model
   double myshadowFactor(GVector& sunpos_eci, GVector& satpos_eci)
  {
           double factor = 1.0;

           double a = 6378.137; //km
           double b = 6356.7523142; //km  6356.7523142

           // for sphere test
          // double a = 6371.0; //km
          // double b = 6371.0; //km  6356.7523142

           double hgt_atm = 50.0; // the hight of atmosphere, this parameter is very important

           //considering atmosphere
           double a_atm = a + hgt_atm; //km
           double b_atm = b*(a_atm/a); //km  6356.7523142

           double Area_earth =0.0, Area_atm = 0.0, Area_solar=0.0;
           int state1 = -1, state2 =-1;

           double x1 =0.0, x2 =0.0, dis0 =0.0, dis1 = 0.0, dis2 = 0.0, thickness;
           double r_sun = 0.0;

           if(hgt_atm == 0.0 )
           {
                state2 = myperspectiveProjection(a,b,sunpos_eci,satpos_eci,r_sun, Area_earth,dis1, dis0);
                Area_solar = 3.14159265357*r_sun*r_sun;

               if(state2 == 0)
               {
                   factor = 1.0;
               }
               else if (state2 == -1)
               {
                   factor =  0.0;
               }
               else
               {
                  factor = Area_earth/Area_solar;
               }

           }
           else  // considering the atmosphere
           {
               state1 = myperspectiveProjection(a_atm,b_atm,sunpos_eci,satpos_eci,r_sun, Area_atm,dis1,dis0);
               state2 = myperspectiveProjection(a,b,sunpos_eci,satpos_eci,r_sun, Area_earth,dis2,dis0);

               double mu1 = 0.0;  // 大气层辐射通过系数, distance = 0
               double mu2 = 1.0;  // distance = 1;

               double a_linear = mu1;
               double b_linear = mu2;

               // so the linear model would be y = 0.6x + 0.3
               double a_exp = mu1;
               double b_exp = std::log(mu2/mu1);

               double a_log = mu1;
               double b_log = (mu2-mu1)/std::log(2.0);

               Area_solar = 3.14159265357*r_sun*r_sun;

               // the thickness of atmosphere
              thickness = dis1 - dis2;

               //线性模型
               //double u1 = (mu2-mu1) * x1 + mu1;
               //double u2 = (mu2-mu1) * x2 + mu1;

               //指数模型
               //double u1 = a_exp*exp(b_exp*x1);
               //double u2 = a_exp*exp(b_exp*x2);

               //对数模型
               //double u1 = a_log + b_log*log(x1+1.0);
               //double u2 = a_log + b_log*log(x2+1.0);

               //double u = a_log + b_log*log(x+1.0);

               // S 曲线模型 Logistic function  y = L/(1 + aexp(bx));
               //double u1 = 1.0/( 1+ exp(-3*x1+1) );
               //double u2 = 1.0/( 1+ exp(-3*x2+1) );

               if(state1 == 0)  // full phase
               {
                   factor = 1.0;
               }
               // partly in the atmosphere
               if((state1 ==1 || state1 == 2 )
                  && state2 == 0
                  )
               {
                   x1 = 1.0;
                   x2 = (thickness - (dis1 - dis0))/thickness;

                   //log
                   //double u1 = a_log + b_log*log(x1+1.0);
                   //double u2 = a_log + b_log*log(x2+1.0);

                   //linear
                   double u1 = (mu2-mu1) * x1 + mu1;
                   double u2 = (mu2-mu1) * x2 + mu1;

                   double coeff = (u1+u2)/2.0;
                   factor = (Area_atm + (Area_solar - Area_atm)*coeff )/Area_solar;
               }

               // totally in the atmosphere
               if( state1 == -1 && state2 == 0 )
               {
                  double u1 =0, u2 =0;
                if(dis2 != 0.0)
                {
                    x1 = (dis0 - dis2)/thickness;
                    x2 = (dis0 - dis2 + 2.0*r_sun )/thickness;

                    //log
                    //double u1 = a_log + b_log*log(x1+1.0);
                    //double u2 = a_log + b_log*log(x2+1.0);

                    //linear
                    u1 = (mu2-mu1) * x1 + mu1;
                    u2 = (mu2-mu1) * x2 + mu1;
                }
                else if (fabs(dis2)< 1.0E-9) // no projection for the solid Earth (outside of the solid Earth), but the satellite is in the atmosphere's projection
                {

                    //thickness = hgt_atm*dis1 / ((a_atm+b_atm)/2.0 );
//                    thickness = 500*r_sun;
//                    x1 = (dis1 - dis0  - 2*r_sun)/thickness;
//                    x2 = (dis1 - dis0)/thickness;
//
//                    //linear
//                    u1 = (mu1-mu2) * x1 + mu2;
//                    u2 = (mu1-mu2) * x2 + mu2;
//
                    u1 = 0.5;
                    u2 = 0.5;

                }

                   factor = (u1+u2)/2.0;
               }
               // partly in the earth and atmosphere
               if( state1 == -1 && (state2 == 1  || state2 == 2) )
               {
                   x1 = (2.0*r_sun - (dis2 - dis0))/thickness;
                   x2 = 0.0;

                   //log
                   //double u1 = a_log + b_log*log(x1+1.0);
                   //double u2 = a_log + b_log*log(x2+1.0);

                   //linear
                   double u1 = (mu2-mu1) * x1 + mu1;
                   double u2 = (mu2-mu1) * x2 + mu1;

                   double coeff = (u1 + u2)/2.0;

                   factor = (Area_earth*coeff )/Area_solar;
               }

               // totally in the earth's shadow, umbra
               // in umbra, the CERES data will capture the atmosphere effect
               if( state2 == -1  )
               {
                   factor = 0.0;
               }

               // the sun is very big, bigger than the atmosphere thickness
               if( (state1 ==1 || state1 == 2)
                  && (state2 ==1 || state2 == 2)
                 )
               {
                   //需要分成3段, 关键确定第2段大气层中的辐射系数
                   double area1 = Area_atm;
                   double area2 = Area_earth - Area_atm;
                   //double area3 = Area_solar - area1 - area2;

                   x1 =  1.0;
                   x2 =  0.0;
                   //对数模型
                   //double u1 = a_log + b_log*log(x1+1.0);
                   //double u2 = a_log + b_log*log(x2+1.0);

                   //linear
                   double u1 = (mu2-mu1) * x1 + mu1;
                   double u2 = (mu2-mu1) * x2 + mu1;

                   factor = (area1*1.0 + area2*(u1+u2)/2.0)/Area_solar;
               }
           }
           return factor;
       }


int main()
{
  // test time: UTC: 2015 1 11 18 33 34
  // shadow function should give: 0.474838
  GVector satpos(-13205.655784525363,21522.519302073124,15446.72240793841); // sat position in km
  GVector sunpos(52727703.80386541,-126017147.89721917,-54630443.258015752); //sun position in km

  double shadow_factor=myshadowFactor(sunpos, satpos);
  printf("shadow factor: %.6f\n",shadow_factor);
  return 0;
}
