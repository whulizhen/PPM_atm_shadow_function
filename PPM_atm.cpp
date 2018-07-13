#include "GVector.hpp"
#include "GMatrix.h"
using namespace gfc;



/*
Modified from https://github.com/dxli/quarticSolver
Author: Dongxu Li, Written for the LibreCAD project: http://librecad.org/
*/
unsigned int quadraticSolver(double * ce,  double * roots)
   //quadratic solver for
   // x^2 + ce[0] x + ce[1] =0
   {
       double discriminant=0.25*ce[0]*ce[0]-ce[1];
       if (discriminant < 0.0 ) return 0;
       roots[0]= -0.5*ce[0] + sqrt(discriminant);
       roots[1]= -ce[0] - roots[0];
       return 2;
   }

   /*
   Modified from https://github.com/dxli/quarticSolver
   Author: Dongxu Li, Written for the LibreCAD project: http://librecad.org/
   */
unsigned int cubicSolver(double * ce, double *roots)
   //cubic equation solver
   // x^3 + ce[0] x^2 + ce[1] x + ce[2] = 0
   {
       // depressed cubic, Tschirnhaus transformation, x= t - b/(3a)
       // t^3 + p t +q =0
       unsigned int ret=0;
       double shift=(1./3)*ce[0];
       double p=ce[1] -shift*ce[0];
       double q=ce[0]*( (2./27)*ce[0]*ce[0]-(1./3)*ce[1])+ce[2];
       //Cardano's method,
       //	t=u+v
       //	u^3 + v^3 + ( 3 uv + p ) (u+v) + q =0
       //	select 3uv + p =0, then,
       //	u^3 + v^3 = -q
       //	u^3 v^3 = - p^3/27
       //	so, u^3 and v^3 are roots of equation,
       //	z^2 + q z - p^3/27 = 0
       //	and u^3,v^3 are,
       //		-q/2 \pm sqrt(q^2/4 + p^3/27)
       //	discriminant= q^2/4 + p^3/27
       //std::cout<<"cubicSolver:: p="<<p<<"\tq="<<q<<std::endl;

       double discriminant= (1./27)*p*p*p+(1./4)*q*q;
       if ( fabs(p)< 1.0e-75)
       {
           ret=1;
           *roots=(q>0)?-pow(q,(1./3)):pow(-q,(1./3));
           *roots -= shift;
           return ret;
       }

       //std::cout<<"cubicSolver:: discriminant="<<discriminant<<std::endl;

       if( discriminant>0.0 )
       {
           double ce2[2]= {q, -1./27*p*p*p},u3[2];
           ret=quadraticSolver(ce2,u3);
           if (! ret )
           { //should not happen
               std::cerr<<"cubicSolver::Error cubicSolver("<<ce[0]<<' '<<ce[1]<<' '<<ce[2]<<")\n";
           }
           ret=1;
           double u,v;
           u= (q<=0) ? pow(u3[0], 1./3): -pow(-u3[1],1./3);
           //u=(q<=0)?pow(-0.5*q+sqrt(discriminant),1./3):-pow(0.5*q+sqrt(discriminant),1./3);
           v=(-1./3)*p/u;

           //std::cout<<"cubicSolver:: u="<<u<<"\tv="<<v<<std::endl;
           //std::cout<<"cubicSolver:: u^3="<<u*u*u<<"\tv^3="<<v*v*v<<std::endl;

           *roots=u+v - shift;
           return ret;
       }
       ret=3;
       std::complex<double> u(q,0),rt[3];
       u=pow(-0.5*u-sqrt(0.25*u*u+p*p*p/27),1./3);
       rt[0]=u-p/(3.*u)-shift;
       std::complex<double> w(-0.5,sqrt(3.)/2);
       rt[1]=u*w-p/(3.*u*w)-shift;
       rt[2]=u/w-p*w/(3.*u)-shift;
       //	std::cout<<"Roots:\n";
       //	std::cout<<rt[0]<<std::endl;
       //	std::cout<<rt[1]<<std::endl;
       //	std::cout<<rt[2]<<std::endl;

       roots[0]=rt[0].real();
       roots[1]=rt[1].real();
       roots[2]=rt[2].real();
       return ret;
   }

   /*
   Modified from https://github.com/dxli/quarticSolver
   Author: Dongxu Li, Written for the LibreCAD project: http://librecad.org/
   */
unsigned int quarticSolver(double * ce, double *roots)
   //quartic solver
   // x^4 + ce[0] x^3 + ce[1] x^2 + ce[2] x + ce[3] = 0
   {
       // x^4 + a x^3 + b x^2 +c x + d = 0
       // depressed quartic, x= t - a/4
       // t^4 + ( b - 3/8 a^2 ) t^2 + (c - a b/2 + a^3/8) t + d - a c /4 + a^2 b/16 - 3 a^4/256 =0
       // t^4 + p t^2 + q t + r =0
       // p= b - (3./8)*a*a;
       // q= c - 0.5*a*b+(1./8)*a*a*a;
       // r= d - 0.25*a*c+(1./16)*a*a*b-(3./256)*a^4
       unsigned int ret=0;
       double shift=0.25*ce[0];
       double shift2=shift*shift;
       double a2=ce[0]*ce[0];
       double p= ce[1] - (3./8)*a2;
       double q= ce[2] + ce[0]*((1./8)*a2 - 0.5*ce[1]);
       double r= ce[3] - shift*ce[2] + (ce[1] - 3.*shift2)*shift2;
       double eps = 1.0E-16;
       //std::cout<<"quarticSolver:: p="<<p<<"\tq="<<q<<"\tr="<<r<<std::endl;

       if (fabs(q) <= eps) // q == 0
       {
           // Biquadratic equations
           double discriminant= 0.25*p*p -r;

           if(discriminant < 0.)
           {
               return 0;
           }

           double t2[2];
           t2[0]=-0.5*p-sqrt(discriminant);  // t2[0] > t2[1]
           t2[1]= -p - t2[0];

           if( fabs(discriminant) < eps )
           {
               roots[0]=sqrt(t2[0])-shift;
               roots[1]= -sqrt(t2[0])-shift;
               return 2;
           }

           if ( t2[0] >= 0. )
           {// four real roots
               roots[0]=sqrt(t2[0])-shift;
               roots[1]= -sqrt(t2[0])-shift;
               roots[2]=sqrt(t2[1])-shift;
               roots[3]= -sqrt(t2[1])-shift;
               return 4;
           }
           if ( t2[1] >= 0.)
           { // two real roots
               roots[0]=sqrt(t2[1])-shift;
               roots[1]= -roots[0]-shift;
               return 2;
           }
           return 0;
       }
       if ( fabs(r)< eps )  // a cubic
       {
           double cubic[3]= {0.,p,q};
           roots[0]=0.;
           ret=1+cubicSolver(cubic,roots+1);
           for(unsigned int i=0; i<ret; i++) roots[i] -= shift;
           return ret;
       }
       double cubic[3]= {2.*p,p*p-4.*r,-q*q},croots[3]={0.0};
       ret = cubicSolver(cubic,croots);
       int rootNum = 0; // this is for roots

       for( int i = 0 ; i< ret ; i++ )
       {
           int num_root = 0;
           double myroots[4] = {0.0};
           double droots[2] = {0.0};
           double ce2[2];

           if( croots[i] < 0.0 || fabs(croots[i])< eps ) // only can be positive
           {
               continue;
           }

           //double u2 = croots[i]; , u can not be ZERO
           double u = sqrt(croots[i]);
           double v = 0.5*(p+croots[i])-0.5*q/u ;
           double w = 0.5*(p+croots[i])+0.5*q/u ;
           ce2[0]=	u;
           ce2[1]= v;
           int ret1=quadraticSolver(ce2,droots);
           for( int j = 0 ; j< ret1 ; j++ )
           {
               if( num_root == 0 )
               {
                   myroots[num_root++] = droots[0];
               }
               else
               {
                   bool exist = false;
                   for(int k = 0 ; k< num_root; k++ )
                   {
                       if( myroots[k] == droots[j])
                       {
                           exist = true;
                           break;
                       }
                   }

                   if(exist == false)
                   {
                       myroots[num_root++] = droots[j];
                   }

               }


           }

           ce2[0]=	-u;
           ce2[1]= w;
           int ret2=quadraticSolver(ce2,droots);
           for( int j = 0 ; j< ret2 ; j++ )
           {
               if( num_root == 0 )
               {
                   myroots[num_root++] = droots[0];
               }
               else
               {
                   bool exist = false;
                   for(int k = 0 ; k< num_root; k++ )
                   {
                       if( myroots[k] == droots[j])
                       {
                           exist = true;
                           break;
                       }
                   }

                   if(exist == false)
                   {
                       myroots[num_root++] = droots[j];
                   }

               }
           }

           //comparing myroots with roots
           for( int j = 0 ; j< num_root ; j++ )
           {
               if( rootNum == 0 )
               {
                   roots[rootNum++] = myroots[0];
               }
               else
               {
                   bool exist = false;
                   for(int k = 0 ; k< rootNum; k++ )
                   {
                       if( roots[k] == myroots[j])
                       {
                           exist = true;
                           break;
                       }
                   }

                   if(exist == false)
                   {
                       roots[rootNum++] = myroots[j];
                   }
               }
           }
       }

       //do the shift
       for(int i = 0 ; i< rootNum; i++ )
       {
           roots[i] -= shift;
       }
       return rootNum;

   }




      double area_hyperbola(double Q1[2], double Q2[2], double Os[2], double Rs,double a, double b, bool in_out, bool x_axis)
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

          //calculate the area of hyperbolic secttion


          double SQ1Q2Oe = 0.5*a*b*fabs( acos(ae[0]/a) - acos(be[0]/a) );
          double S1 = SQ1Q2Oe - TQ1Q2Oe;

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
          double s1 = b*( Q1[0]*sqrt(  (Q1[0]*Q1[0]/a/a) -1.0 ) - a*acosh(Q1[0]/a) );
          double s2 = b*( Q2[0]*sqrt(  (Q2[0]*Q2[0]/a/a) -1.0 ) - a*acosh(Q2[0]/a) );
          double ss = s1<=s2?s1:s2;
          double sl = s1>=s2?s1:s2;
          double S2 =0.0;
          if( Q1[1]*Q2[1] < 0.0 )  //2个交点在x轴不同侧
          {
              double k = ( Q2[1] - Q1[1] )/( Q2[0] - Q1[0]);
              double m = Q1[1] - k *Q1[0];
              double x_m = - m/k;

              S2 = ss + (sl - ss)/2.0 - 0.5*fabs((Q2[0]- x_m)*Q2[1]) + 0.5*fabs((x_m - Q1[0])*Q1[1]);
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

          double SQ1Q2Oe = 0.5*a*b*fabs( acos(ae[0]/a) - acos(be[0]/a) );

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

       GVector o, u, v, n;
       double f = 1000.0;
       n = r - rs;
       double dis_sat_earth = r.norm();
       double dis_sat_sun = n.norm();
       n.normalise();

       //A: first test if the satellite is in the front of earth,
       // if it is in the front of earth, it is always full phase
       // otherwise,using the photogrammetry method
       // A= diag{1/a2,1/a2, 1/b2 }, A^{-1} = diag{a2, a2, b2}
       double nAn = n.x*n.x*a2 + n.y*n.y*a2 + + n.z*n.z*b2;
       double s1 = sqrt(1.0/nAn);
       double s2 = - s1;
       double rtn = dotproduct(r, n);
       double ntn = dotproduct(n, n);
       double t1 = (rtn - 1.0/s1)/ntn;
       double t2 = (rtn - 1.0/s2)/ntn;

       //GVector xp1(s1*n.x/a2, s1*n.y/b2, s1*n.z/b2 );
       //GVector xp2(s2*n.x/a2, s2*n.y/b2, s2*n.z/b2 );

       GVector xs1 = r - t1*n;
       GVector xs2 = r - t2*n;

       double ds1 = (xs1-rs).norm();
       double ds2 = (xs2-rs).norm();

       double t =0.0;
       t = (ds1 <= ds2) ? ds1:ds2;
       ds2 = (ds1 >= ds2)? ds1:ds2;
       ds1 = t;

       if(dis_sat_sun < ds1)  // full phase
       {
           state = 0;
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
       //u = r - rtn*n;
       //u = n - 1.0/rtn*r;
       u = rtn*n - r;
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

       //C: calculate M
       GMatrix M(3,3);//A(3,3);
       t = (r.x*r.x/a2 + r.y*r.y/a2 + r.z*r.z/b2 - 1.0);
       M(0,0) = r.x*r.x/a2/a2 - t/a2;
       M(0,1) = r.x*r.y/a2/a2;
       M(0,2) = r.x*r.z/a2/b2;
       M(1,0) = M(0,1);
       M(1,1) = r.y*r.y/a2/a2 - t/a2;
       M(1,2) = r.y*r.z/a2/b2;
       M(2,0) = M(0,2);
       M(2,1) = M(1,2);
       M(2,2) = r.z*r.z/b2/b2 - t/b2;

       //D: calculate K
       double K[6] = {0.0};
       GMatrix tu(3,1), tv(3,1),tn(3,1), to(3,1),tr(3,1);
       tu[0] = u.x;tu[1] = u.y;tu[2] = u.z;
       tv[0] = v.x;tv[1] = v.y;tv[2] = v.z;
       //to[0] = o.x;to[1] = o.y;to[2] = o.z;
       //tr[0] = r.x;tr[1] = r.y;tr[2] = r.z;
       tn[0] = n.x;tn[1] = n.y;tn[2] = n.z;

       ((~tu)*M*tu).getData(K); // K[0]
       ((~tu)*M*tv).getData(K+1); // K[1]
       ((~tv)*M*tv).getData(K+2); // K[2]
       (-2.0*f*(~tu*M*tn)).getData(K+3); //K[3]
       (-2.0*f*(~tv*M*tn)).getData(K+4);  //K[4]
       (f*f*(~tn*M*tn)).getData(K+5); //K[5]




       // just make the numbers larger for visualisation
       for( int i =0; i< 6; i++ )
       {
           K[i] = K[i]*1.0E6;
       }

//        double CEOF[5]={K[0]/K[5], 2*K[1]/K[5], K[2]/K[5], K[3]/K[5],K[4]/K[5] };
//        double XC[2] = { (CEOF[1]*CEOF[4] - 2*CEOF[2]*CEOF[3])/(4*CEOF[0]*CEOF[2] - CEOF[1]*CEOF[1]), (CEOF[1]*CEOF[3] - 2*CEOF[0]*CEOF[4])/(4*CEOF[0]*CEOF[2] - CEOF[1]*CEOF[1])    };


       //E: calculate the eigen value and eigen vector of matrix B
       t =  sqrt( (K[0]+K[2])*(K[0]+K[2]) - 4.0*(K[0]*K[2] - K[1]*K[1]) ) ;
       double lambda1 = (K[0]+K[2]+t)/2.0;
       double lambda2 = (K[0]+K[2]-t)/2.0;

       //get the eigen vector of B, i.e. matrix Q
       double  r1[2]={0.0,1.0},r2[2]={1.0,0.0};
       bool sol1 = false, sol2 =false;
       if( fabs(lambda1 - K[0])<1.0E-12)
       {
           r1[1] = 0.0;
           if( (lambda1 - K[0])*K[1] <0.0 ) // 异号
           {
               r1[0] = -1.0;
           }
           else
           {
               r1[0] = 1.0;
           }
           sol1 = true;
       }

       if( fabs(lambda2 - K[2])<1.0E-12)
       {
           r2[0] = 0.0;
           if( (lambda2 - K[2])*K[1] <0.0 ) // 异号
           {
               r2[1] = -1.0;
           }
           else
           {
               r2[1] = 1.0;
           }
           sol2 = true;
       }

       if(sol1 ==false)
       {
           r1[0] = K[1]/(lambda1 - K[0]);  // r1 for lambda1
       }

       if(sol2 ==false)
       {
           r2[1] = K[1]/(lambda2 - K[2]);  // r2 for lambda2
       }

       //get the unit vector
       t = sqrt(r1[0]*r1[0]+r1[1]*r1[1]);
       r1[0] = r1[0]/t;r1[1] = r1[1]/t;

       t= sqrt(r2[0]*r2[0]+r2[1]*r2[1]);
       r2[0] = r2[0]/t;r2[1] = r2[1]/t;

       //the larger eigen value is lambda1
       if( fabs(lambda1) <= fabs(lambda2) )  // swap lamda1 and lamda2, with r1 and r2
       {
           t = lambda1;
           lambda1 = lambda2;
           lambda2 = t;
           double r[2] ={r1[0],r1[1]};
           r1[0] = r2[0]; r1[1] = r2[1];
           r2[0] = r[0]; r2[1] = r[1];
       }

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
       num_of_solution = quarticSolver(ce,XX);

       if( num_of_solution == 3 || num_of_solution == 4)
       {
           printf("WARNING: shadowfactor, %d intersections!\n", num_of_solution);
           exit(0);
       }

       // calculate the coordinates of the intersections betweent the circle and the conical curve in transformed frame
       // in the new frame, the origin is at the centre of the projection of the Earth
       double Q1[2]={ (1.0 - XX[0]*XX[0])/(1.0 + XX[0]*XX[0])*R0 + 0.5*tx , 2*XX[0]/(1.0 + XX[0]*XX[0])*R0 + 0.5*ty};
       double Q2[2]={ (1.0 - XX[1]*XX[1])/(1.0 + XX[1]*XX[1])*R0 + 0.5*tx , 2*XX[1]/(1.0 + XX[1]*XX[1])*R0 + 0.5*ty};

       //calculate the sun's projection centre Os (PSC) in the transformed frame
       double Os[2] = {0.5*tx , 0.5*ty };
       //double Os[2] = {0.0,0.0};
       double PEC_new[2] = {0.0,0.0};
       // the projection of the earth's centre Pe in the transformed frame is obtained by converting PEC_isf into the rotated and translated frame
       PEC_new[0] = PEC_isf.x*r1[0] + r1[1]*PEC_isf.y + 0.5*tx;
       PEC_new[1] = PEC_isf.x*r2[0] + r2[1]*PEC_isf.y + 0.5*ty;

       //PEC_new[0] = (PEC_isf.x-0.5*tx)*r1[0] + (PEC_isf.y-0.5*ty)*r2[0];
       //PEC_new[1] = (PEC_isf.x-0.5*tx)*r1[1] + (PEC_isf.y-0.5*ty)*r2[1];


       //double boundary_intersection[2] ={0.0};
       //PEC_new[0] = 0.0;
       //PEC_new[1] = 0.0;

       //double mytest = (lambda1*Os[0]*Os[0] + lambda2*Os[1]*Os[1] )/(OM-K[5]);

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

           area_bright = area_hyperbola( Q1, Q2, Os, R0, a, b, in_out, x_axis);
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

           double EC_intersection1[2]={0.0},EC_intersection2[2]={0.0}, EC[2]={0.0};
           double Boundary_earth[2] = {0.0}, Boundary_atm[2] = {0.0};
           double PEC_new[2] ={0.0};
           double sun_intersection[2]={0.0};

           double x1 =0.0, x2 =0.0, dis0 =0.0, dis1 = 0.0, dis2 = 0.0, thickness;
           double r_sun = 0.0;

           if(hgt_atm == 0.0 )
           {
               //state2 = perspectiveProjection(a,b,sunpos_eci,satpos_eci,r_sun, Area_earth,EC_intersection2, EC);
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
               state2 = myperspectiveProjection(a,b,sunpos_eci,satpos_eci,r_sun, Area_earth,dis2,dis0);
               state1 = myperspectiveProjection(a_atm,b_atm,sunpos_eci,satpos_eci,r_sun, Area_atm,dis1,dis0);

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
                   x1 = (dis0 - dis2)/thickness;
                   x2 = (dis0 - dis2 + 2.0*r_sun )/thickness;

                   //log
                   //double u1 = a_log + b_log*log(x1+1.0);
                   //double u2 = a_log + b_log*log(x2+1.0);

                   //linear
                   double u1 = (mu2-mu1) * x1 + mu1;
                   double u2 = (mu2-mu1) * x2 + mu1;

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
