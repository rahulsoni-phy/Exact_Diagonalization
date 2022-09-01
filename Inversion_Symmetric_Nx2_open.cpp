#include <iostream>
#include <math.h>
#include <time.h>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include "tensor.h"
#include <assert.h>
#include <algorithm>
#define PI 3.14159265359

int main(int argc, char** argv){

int Lx, Ly, Mat_Size, Mat_Half_Size;
double hop, c_pot, l_coup, B_mag, B1_mag;

string connection_out;

//----------------------Inputs-------------------//
Lx=atoi(argv[1]);
Ly=2;

l_coup=atof(argv[2]);
B_mag=atof(argv[3]);
B1_mag=atof(argv[4]);
c_pot=atof(argv[5]);

connection_out=argv[6];
ofstream connec_out(connection_out.c_str());

//-----------------------------------------------//
int r, r1, r2, r3, r4, r5, r6;
Mat_Size=6*Lx*Ly-8;
Mat_Half_Size=3*Lx*Ly-4;
hop=1.0;
cout<<"This code generates OBCxOBC inversion symmetric long-range Dice connections"<<endl;
cout<<"Total number of sites in the lattice ="<<Mat_Half_Size<<endl;
cout<<"Rashba-SOC strength ="<<l_coup<<endl;
cout<<"Magnetic field strength on site type-1="<<B_mag<<endl;
cout<<"Magnetic field strength on site type-2="<<B1_mag<<endl;
cout<<"Onsite-energy="<<c_pot<<endl;

Mat_2_Complex_doub C_mat;
C_mat.resize(Mat_Size);

for(int i=0;i<Mat_Size;i++){
    C_mat[i].resize(Mat_Size);
    for(int j=0;j<Mat_Size;j++){
        C_mat[i][j].real(0.0);
        C_mat[i][j].imag(0.0);
    }
}

//C_mat[i][j] is C_mat[Create at i][Annihilate at j]

int MHS = Mat_Half_Size;
//---------------------------------Adding Kinetic Term----------------------------//
for(int sp=0;sp<2;sp++){
for(int ap=0;ap<3;ap++){

		if(ap==1){
			C_mat[r][r1].real( (-1.0)*hop );
			C_mat[r][r1].imag(0.0);
			C_mat[r1][r].real( (-1.0)*hop );
			C_mat[r1][r].imag(0.0);

			C_mat[r2][r].real( (-1.0)*hop );
                        C_mat[r2][r].imag(0.0);
                        C_mat[r][r2].real( (-1.0)*hop );
                        C_mat[r][r2].imag(0.0);

			if(iy!=Ly-1){
				C_mat[r][r3].real( (-1.0)*hop );
	                        C_mat[r][r3].imag(0.0);
        	                C_mat[r3][r].real( (-1.0)*hop );
                	        C_mat[r3][r].imag(0.0);
			}

			if(ix!=Lx-1){
				C_mat[r][r4].real( (-1.0)*hop );
        	                C_mat[r][r4].imag(0.0);
                	        C_mat[r4][r].real( (-1.0)*hop );
                        	C_mat[r4][r].imag(0.0);
			}

			if(ix!=0){
                                C_mat[r][r5].real( (-1.0)*hop );
                                C_mat[r][r5].imag(0.0);
                                C_mat[r5][r].real( (-1.0)*hop );
                                C_mat[r5][r].imag(0.0);
                        }

                        if(iy!=0){
                                C_mat[r][r6].real( (-1.0)*hop );
                                C_mat[r][r6].imag(0.0);
                                C_mat[r6][r].real( (-1.0)*hop );
                                C_mat[r6][r].imag(0.0);
                        }
		}
}
}

//--------------------------Adding Chemical Potential Term and Magnetic Field Term-------------------------//
for(int sp=0;sp<2;sp++){
for(int ap=0;ap<3;ap++){
        for(int ix=0;ix<Lx;ix++){
        for(int iy=0;iy<Ly;iy++){
		r  = sp*MHS + (ap   + 3*ix     + 3*Lx*iy );

		if(sp==0){
		if(ap==1){
		C_mat[r][r].real( (-1.0)*B1_mag - (1.0)*c_pot );
		C_mat[r][r].imag(0.0);	
		}
		else{
		C_mat[r][r].real( (-1.0)*B_mag);
                C_mat[r][r].imag(0.0);
		}
	
	}

		if(sp==1){
		if(ap==1){
                C_mat[r][r].real( (1.0)*B1_mag - (1.0)*c_pot );
                C_mat[r][r].imag(0.0);
                }else{
        	C_mat[r][r].real( (1.0)*B_mag);
        	C_mat[r][r].imag(0.0);        
        	}
		}
	}
	}

}
}

//------------------------------------Adding Rashba-SOC--------------------------------------------------//
for(int ap=0;ap<3;ap++){
        for(int ix=0;ix<Lx;ix++){
        for(int iy=0;iy<Ly;iy++){
                        r  = (ap   + 3*ix     + 3*Lx*iy );
                        r1 = (ap-1 + 3*ix     + 3*Lx*iy );
                        r2 = (ap+1 + 3*ix     + 3*Lx*iy );
			r3 = (ap-1 + 3*(ix+1) + 3*Lx*iy );
			r4 = (ap-1 + 3*ix     + 3*Lx*(iy+1) );
			r5 = (ap+1 + 3*(ix-1) + 3*Lx*iy );
			r6 = (ap+1 + 3*ix     + 3*Lx*(iy-1) );

		if(ap==1){
			C_mat[r1][r+MHS].real( 0.0 ); 			//h_{i1,j}_{dn->up}  = -i*lambda (c^{d}_{i1,up}c_{j,dn} - c^{d}_{j,dn}c_{i1,up})
			C_mat[r1][r+MHS].imag( (-1.0)*l_coup );
			C_mat[r+MHS][r1].real( 0.0 );
                        C_mat[r+MHS][r1].imag( (1.0)*l_coup );

			C_mat[r1+MHS][r].real( 0.0 );			//h_{i1,j}_{up->dn}  = -i*lambda (c^{d}_{i1,dn}c_{j,up} - c^{d}_{j,up}c_{i1,dn})
                        C_mat[r1+MHS][r].imag( (-1.0)*l_coup );
                        C_mat[r][r1+MHS].real( 0.0 );
                        C_mat[r][r1+MHS].imag( (1.0)*l_coup );

			C_mat[r2][r+MHS].real( 0.0 );                   //h_{i4,j}_{dn->up}  = -i*lambda (c^{d}_{i4,up}c_{j,dn} - c^{d}_{j,dn}c_{i4,up})
                        C_mat[r2][r+MHS].imag( (-1.0)*l_coup );
                        C_mat[r+MHS][r2].real( 0.0 );
                        C_mat[r+MHS][r2].imag( (1.0)*l_coup );

                        C_mat[r2+MHS][r].real( 0.0 );                   //h_{i4,j}_{up->dn}  = -i*lambda (c^{d}_{i4,dn}c_{j,up} - c^{d}_{j,up}c_{i4,dn})
                        C_mat[r2+MHS][r].imag( (-1.0)*l_coup );
                        C_mat[r][r2+MHS].real( 0.0 );
                        C_mat[r][r2+MHS].imag( (1.0)*l_coup );


			if(ix!=Lx-1){
			C_mat[r3][r+MHS].real( (-1.0)*l_coup*sin( (2.0*PI)/3.0 ) );       //h_{i3,j}_{dn->up}  = -i*lambda (e^{-i*2pi/3}c^{d}_{i3,up}c_{j,dn} - e^{i*2pi/3}c^{d}_{j,dn}c_{i3,up})
                        C_mat[r3][r+MHS].imag( (-1.0)*l_coup*cos( (2.0*PI)/3.0 ) );
                        C_mat[r+MHS][r3].real( (-1.0)*l_coup*sin( (2.0*PI)/3.0 ) );
                        C_mat[r+MHS][r3].imag( (1.0)*l_coup*cos( (2.0*PI)/3.0 ) );

                        C_mat[r3+MHS][r].real( (1.0)*l_coup*sin( (2.0*PI)/3.0 ) );        //h_{i3,j}_{up->dn}  = -i*lambda (e^{i*2pi/3}c^{d}_{i3,dn}c_{j,up} - e^{-i*2pi/3}c^{d}_{j,up}c_{i3,dn})
                        C_mat[r3+MHS][r].imag( (-1.0)*l_coup*cos( (2.0*PI)/3.0 ) );
                        C_mat[r][r3+MHS].real( (1.0)*l_coup*sin( (2.0*PI)/3.0 ) );
                        C_mat[r][r3+MHS].imag( (1.0)*l_coup*cos( (2.0*PI)/3.0 ) );
			}


			if(iy!=Ly-1){
			C_mat[r4][r+MHS].real( (1.0)*l_coup*sin( (2.0*PI)/3.0 ) );        //h_{i5,j}_{dn->up}  = -i*lambda (e^{i*2pi/3}c^{d}_{i5,up}c_{j,dn} - e^{-i*2pi/3}c^{d}_{j,dn}c_{i5,up})
                        C_mat[r4][r+MHS].imag( (-1.0)*l_coup*cos( (2.0*PI)/3.0 ) );
                        C_mat[r+MHS][r4].real( (1.0)*l_coup*sin( (2.0*PI)/3.0 ) );
                        C_mat[r+MHS][r4].imag( (1.0)*l_coup*cos( (2.0*PI)/3.0 ) );

                        C_mat[r4+MHS][r].real( (-1.0)*l_coup*sin( (2.0*PI)/3.0 ) );       //h_{i5,j}_{up->dn}  = -i*lambda (e^{-i*2pi/3}c^{d}_{i5,dn}c_{j,up} - e^{i*2pi/3}c^{d}_{j,up}c_{i5,dn})
                        C_mat[r4+MHS][r].imag( (-1.0)*l_coup*cos( (2.0*PI)/3.0 ) );
                        C_mat[r][r4+MHS].real( (-1.0)*l_coup*sin( (2.0*PI)/3.0 ) );
                        C_mat[r][r4+MHS].imag( (1.0)*l_coup*cos( (2.0*PI)/3.0 ) );
			}

			if(ix!=0){
			C_mat[r5][r+MHS].real( (-1.0)*l_coup*sin( (2.0*PI)/3.0 ) );       //h_{i6,j}_{dn->up}  = -i*lambda (e^{-i*2pi/3}c^{d}_{i6,up}c_{j,dn} - e^{i*2pi/3}c^{d}_{j,dn}c_{i6,up})
                        C_mat[r5][r+MHS].imag( (-1.0)*l_coup*cos( (2.0*PI)/3.0 ) );
                        C_mat[r+MHS][r5].real( (-1.0)*l_coup*sin( (2.0*PI)/3.0 ) );
                        C_mat[r+MHS][r5].imag( (1.0)*l_coup*cos( (2.0*PI)/3.0 ) );

                        C_mat[r5+MHS][r].real( (1.0)*l_coup*sin( (2.0*PI)/3.0 ) );        //h_{i6,j}_{up->dn}  = -i*lambda (e^{i*2pi/3}c^{d}_{i6,dn}c_{j,up} - e^{-i*2pi/3}c^{d}_{j,up}c_{i6,dn})
                        C_mat[r5+MHS][r].imag( (-1.0)*l_coup*cos( (2.0*PI)/3.0 ) );
                        C_mat[r][r5+MHS].real( (1.0)*l_coup*sin( (2.0*PI)/3.0 ) );
                        C_mat[r][r5+MHS].imag( (1.0)*l_coup*cos( (2.0*PI)/3.0 ) );
			}

			if(iy!=0){
			C_mat[r6][r+MHS].real( (1.0)*l_coup*sin( (2.0*PI)/3.0 ) );        //h_{i2,j}_{dn->up}  = -i*lambda (e^{i*2pi/3}c^{d}_{i2,up}c_{j,dn} - e^{-i*2pi/3}c^{d}_{j,dn}c_{i2,up})
                        C_mat[r6][r+MHS].imag( (-1.0)*l_coup*cos( (2.0*PI)/3.0 ) );
                        C_mat[r+MHS][r6].real( (1.0)*l_coup*sin( (2.0*PI)/3.0 ) );
                        C_mat[r+MHS][r6].imag( (1.0)*l_coup*cos( (2.0*PI)/3.0 ) );

                        C_mat[r6+MHS][r].real( (-1.0)*l_coup*sin( (2.0*PI)/3.0 ) );       //h_{i2,j}_{up->dn}  = -i*lambda (e^{-i*2pi/3}c^{d}_{i2,dn}c_{j,up} - e^{i*2pi/3}c^{d}_{j,up}c_{i2,dn})
                        C_mat[r6+MHS][r].imag( (-1.0)*l_coup*cos( (2.0*PI)/3.0 ) );
                        C_mat[r][r6+MHS].real( (-1.0)*l_coup*sin( (2.0*PI)/3.0 ) );
                        C_mat[r][r6+MHS].imag( (1.0)*l_coup*cos( (2.0*PI)/3.0 ) );
			}
		}
	

	}
	}
}


/*
for(int i=0;i<L;i++){
for(int j=0;j<L;j++){
cout<<C_mat[i][j]<<endl;
}
}
*/


for (int i=0;i<Mat_Size;i++){
    for (int j=0;j<Mat_Size;j++)
    {
    connec_out<<C_mat[i][j]<<" ";
    }
    connec_out<<endl;
}


double sign_B1, sign_B2;

sign_B1=B_mag/abs(B_mag);
sign_B2=B1_mag/abs(B1_mag);

connec_out<<sign_B1<<"	"<<sign_B2<<"	"<<sign_B1<<endl;

return 0;
}


