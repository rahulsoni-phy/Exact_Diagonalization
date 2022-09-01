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
double hop, hop_pbc, c_pot, l_coup, B1_mag, B2_mag;

string connection_out;

//----------------------Inputs-------------------//
Lx=atoi(argv[1]);
Ly=atoi(argv[2]);

l_coup=atof(argv[3]);

B1_mag=atof(argv[4]);

B2_mag=atof(argv[5]);

c_pot=0.0;

connection_out=argv[6];
ofstream connec_out(connection_out.c_str());


bool Periodic=false;
bool Half_Periodic=false;

//-----------------------------------------------//
int r, r1, r2, r3, r4, r5, r6;
Mat_Size=6*Lx*Ly;
Mat_Half_Size=3*Lx*Ly;
hop=1.0;
cout<<"Total number of unit cells in the lattice ="<<Lx*Ly<<endl;
cout<<"Total number of atoms in the lattice ="<<3*Lx*Ly<<endl;
cout<<"Rashba-SOC strength ="<<l_coup<<endl;
cout<<"Magnetic field strength on six co-ordination site="<<B2_mag<<endl;
cout<<"Magnetic field strength on three co-ordination site="<<B1_mag<<endl;

Mat_2_Complex_doub C_mat;
C_mat.resize(Mat_Size);

for(int i=0;i<Mat_Size;i++){
    C_mat[i].resize(Mat_Size);
    for(int j=0;j<Mat_Size;j++){
        C_mat[i][j].real(0.0);
        C_mat[i][j].imag(0.0);
    }
}

//---------------------------------Adding Kinetic Term----------------------------//
for(int sp=0;sp<2;sp++){
for(int ap=0;ap<3;ap++){
	for(int ix=0;ix<Lx;ix++){
	for(int iy=0;iy<Ly;iy++){
			r  = sp*Mat_Half_Size + (ap   + 3*ix     + 3*Lx*iy );
                        r1 = sp*Mat_Half_Size + (ap-1 + 3*ix     + 3*Lx*iy );
                        r2 = sp*Mat_Half_Size + (ap-1 + 3*(ix+1) + 3*Lx*iy );
                        r3 = sp*Mat_Half_Size + (ap-1 + 3*ix     + 3*Lx*(iy+1) );

		if(ap!=0){
		if(ix!=Lx-1){
		if(iy!=Ly-1){
			C_mat[r][r1].real( (-1.0)*hop );
			C_mat[r][r1].imag(0.0);
			C_mat[r1][r].real( (-1.0)*hop );
			C_mat[r1][r].imag(0.0);

			
			C_mat[r2][r].real( (-1.0)*hop );
			C_mat[r2][r].imag(0.0);
			C_mat[r][r2].real( (-1.0)*hop );
			C_mat[r][r2].imag(0.0);

			C_mat[r3][r].real( (-1.0)*hop );
                        C_mat[r3][r].imag(0.0);
                        C_mat[r][r3].real( (-1.0)*hop );
                        C_mat[r][r3].imag(0.0);
                }
                }
		}

		if(ap!=0){
                if(ix==Lx-1){
                if(iy!=Ly-1){
                        C_mat[r][r1].real( (-1.0)*hop );
                        C_mat[r][r1].imag(0.0);
                        C_mat[r1][r].real( (-1.0)*hop );
                        C_mat[r1][r].imag(0.0);

                        C_mat[r3][r].real( (-1.0)*hop );
                        C_mat[r3][r].imag(0.0);
                        C_mat[r][r3].real( (-1.0)*hop );
                        C_mat[r][r3].imag(0.0);
                }
                }
                }

		if(ap!=0){
                if(ix!=Lx-1){
                if(iy==Ly-1){
                        C_mat[r][r1].real( (-1.0)*hop );
                        C_mat[r][r1].imag(0.0);
                        C_mat[r1][r].real( (-1.0)*hop );
                        C_mat[r1][r].imag(0.0);

                        C_mat[r2][r].real( (-1.0)*hop );
                        C_mat[r2][r].imag(0.0);
                        C_mat[r][r2].real( (-1.0)*hop );
                        C_mat[r][r2].imag(0.0);
                }
                }
                }

		if(ap!=0){
                if(ix==Lx-1){
                if(iy==Ly-1){
                        C_mat[r][r1].real( (-1.0)*hop );
                        C_mat[r][r1].imag(0.0);
                        C_mat[r1][r].real( (-1.0)*hop );
                        C_mat[r1][r].imag(0.0);

                }
                }
                }
		
	}
	}
}
}


if(Periodic==true){
for(int sp=0;sp<2;sp++){
for(int ap=0;ap<3;ap++){
        for(int ix=0;ix<Lx;ix++){
        for(int iy=0;iy<Ly;iy++){

		r  = sp*Mat_Half_Size + (ap   + 3*ix     + 3*Lx*iy );

		if(Half_Periodic==false){
                if(iy==Ly-1){
                  if(ap!=0){
                        r3 = sp*Mat_Half_Size + (ap-1 + 3*ix );

			C_mat[r3][r].real( (-1.0)*hop );
                        C_mat[r3][r].imag(0.0);
                        C_mat[r][r3].real( (-1.0)*hop );
                        C_mat[r][r3].imag(0.0);
                  }
		}
		}
		
		if(ix==Lx-1){
                  if(ap!=0){
                        r2 = sp*Mat_Half_Size + (ap-1 + 3*Lx*iy );
                        
                        C_mat[r2][r].real( (-1.0)*hop );
                        C_mat[r2][r].imag(0.0);
                        C_mat[r][r2].real( (-1.0)*hop );
                        C_mat[r][r2].imag(0.0);
                  }
                }

	}
	}
}
}
}

//--------------------------Adding Chemical Potential Term and Magnetic Field Term-------------------------//
for(int sp=0;sp<2;sp++){
for(int ap=0;ap<3;ap++){
        for(int ix=0;ix<Lx;ix++){
        for(int iy=0;iy<Ly;iy++){
		r  = sp*Mat_Half_Size + (ap   + 3*ix     + 3*Lx*iy );

		if(sp==0){
		if(ap==1){
		C_mat[r][r].real( (-1.0)*B2_mag - (1.0)*c_pot );
		C_mat[r][r].imag(0.0);	
		}
		else{
		C_mat[r][r].real( (-1.0)*B1_mag);
                C_mat[r][r].imag(0.0);
		}
		}

		if(sp==1){
		if(ap==1){
                C_mat[r][r].real( (1.0)*B2_mag - (1.0)*c_pot );
                C_mat[r][r].imag(0.0);
                }else{
        	C_mat[r][r].real( (1.0)*B1_mag);
        	C_mat[r][r].imag(0.0);        
        	}
		}
	}
	}

}
}

//------------------------------------Adding Rashba-SOC--------------------------------------------------//
int MHS = Mat_Half_Size;

for(int ap=0;ap<3;ap++){
        for(int ix=0;ix<Lx;ix++){
        for(int iy=0;iy<Ly;iy++){
                        r  = (ap   + 3*ix     + 3*Lx*iy );
                        r1 = (ap-1 + 3*ix     + 3*Lx*iy );
                        r2 = (ap-1 + 3*(ix+1) + 3*Lx*iy );
                        r3 = (ap-1 + 3*ix     + 3*Lx*(iy+1) );

		if(ap==1){
		if(ix!=Lx-1){
		if(iy!=Ly-1){
			C_mat[r1+MHS][r].real( (-1.0)*l_coup*sin(PI/3) );
                        C_mat[r1+MHS][r].imag( (-1.0)*l_coup*cos(PI/3) );
                        C_mat[r][r1+MHS].real( (-1.0)*l_coup*sin(PI/3) );
                        C_mat[r][r1+MHS].imag( (1.0)*l_coup*cos(PI/3) );

			C_mat[r1][r+MHS].real( (1.0)*l_coup*sin(PI/3)  );
                        C_mat[r1][r+MHS].imag( (-1.0)*l_coup*cos(PI/3) );
                        C_mat[r+MHS][r1].real( (1.0)*l_coup*sin(PI/3) );
                        C_mat[r+MHS][r1].imag( (1.0)*l_coup*cos(PI/3) );

			C_mat[r2+MHS][r].real( (1.0)*l_coup*sin(PI/3) );
                        C_mat[r2+MHS][r].imag( (-1.0)*l_coup*cos(PI/3) );
                        C_mat[r][r2+MHS].real( (1.0)*l_coup*sin(PI/3) );
                        C_mat[r][r2+MHS].imag( (1.0)*l_coup*cos(PI/3) );
                 
                        C_mat[r2][r+MHS].real( (-1.0)*l_coup*sin(PI/3)  );
                        C_mat[r2][r+MHS].imag( (-1.0)*l_coup*cos(PI/3) );
                        C_mat[r+MHS][r2].real( (-1.0)*l_coup*sin(PI/3) );
                        C_mat[r+MHS][r2].imag( (1.0)*l_coup*cos(PI/3) );
			
			C_mat[r3+MHS][r].real( 0.0 );
                        C_mat[r3+MHS][r].imag( (1.0)*l_coup );
                        C_mat[r][r3+MHS].real( 0.0 );
                        C_mat[r][r3+MHS].imag( (-1.0)*l_coup );
                 
                        C_mat[r3][r+MHS].real( 0.0 );
                        C_mat[r3][r+MHS].imag( (1.0)*l_coup );
                        C_mat[r+MHS][r3].real( 0.0 );
                        C_mat[r+MHS][r3].imag( (-1.0)*l_coup );
		}
		}
		}

		if(ap==2){
                if(ix!=Lx-1){
                if(iy!=Ly-1){
                        C_mat[r+MHS][r1].real( (-1.0)*l_coup*sin(PI/3) );
                        C_mat[r+MHS][r1].imag( (-1.0)*l_coup*cos(PI/3) );
                        C_mat[r1][r+MHS].real( (-1.0)*l_coup*sin(PI/3) );
                        C_mat[r1][r+MHS].imag( (1.0)*l_coup*cos(PI/3) );
                 
                        C_mat[r][r1+MHS].real( (1.0)*l_coup*sin(PI/3)  );
                        C_mat[r][r1+MHS].imag( (-1.0)*l_coup*cos(PI/3) );
                        C_mat[r1+MHS][r].real( (1.0)*l_coup*sin(PI/3) );
                        C_mat[r1+MHS][r].imag( (1.0)*l_coup*cos(PI/3) );

                        C_mat[r+MHS][r2].real( (1.0)*l_coup*sin(PI/3) );
                        C_mat[r+MHS][r2].imag( (-1.0)*l_coup*cos(PI/3) );
                        C_mat[r2][r+MHS].real( (1.0)*l_coup*sin(PI/3) );
                        C_mat[r2][r+MHS].imag( (1.0)*l_coup*cos(PI/3) );

                        C_mat[r][r2+MHS].real( (-1.0)*l_coup*sin(PI/3)  );
                        C_mat[r][r2+MHS].imag( (-1.0)*l_coup*cos(PI/3) );
                        C_mat[r2+MHS][r].real( (-1.0)*l_coup*sin(PI/3) );
                        C_mat[r2+MHS][r].imag( (1.0)*l_coup*cos(PI/3) );

                        C_mat[r+MHS][r3].real( 0.0 );
                        C_mat[r+MHS][r3].imag( (1.0)*l_coup );
                        C_mat[r][r+MHS].real( 0.0 );
                        C_mat[r][r+MHS].imag( (-1.0)*l_coup );

                        C_mat[r][r3+MHS].real( 0.0 );
                        C_mat[r][r3+MHS].imag( (1.0)*l_coup );
                        C_mat[r3+MHS][r].real( 0.0 );
                        C_mat[r3+MHS][r].imag( (-1.0)*l_coup );
                }
                }
                }

		if(ap==1){
                if(ix==Lx-1){
                if(iy!=Ly-1){
	                C_mat[r1+MHS][r].real( (-1.0)*l_coup*sin(PI/3) );
                        C_mat[r1+MHS][r].imag( (-1.0)*l_coup*cos(PI/3) );
                        C_mat[r][r1+MHS].real( (-1.0)*l_coup*sin(PI/3) );
                        C_mat[r][r1+MHS].imag( (1.0)*l_coup*cos(PI/3) );

                        C_mat[r1][r+MHS].real( (1.0)*l_coup*sin(PI/3)  );
                        C_mat[r1][r+MHS].imag( (-1.0)*l_coup*cos(PI/3) );
                        C_mat[r+MHS][r1].real( (1.0)*l_coup*sin(PI/3) );
                        C_mat[r+MHS][r1].imag( (1.0)*l_coup*cos(PI/3) );

			C_mat[r3+MHS][r].real( 0.0 );
                        C_mat[r3+MHS][r].imag( (1.0)*l_coup );
                        C_mat[r][r3+MHS].real( 0.0 );
                        C_mat[r][r3+MHS].imag( (-1.0)*l_coup );

                        C_mat[r3][r+MHS].real( 0.0 );
                        C_mat[r3][r+MHS].imag( (1.0)*l_coup );
                        C_mat[r+MHS][r3].real( 0.0 );
                        C_mat[r+MHS][r3].imag( (-1.0)*l_coup );
		}
                }
                }

		if(ap==2){
                if(ix==Lx-1){
                if(iy!=Ly-1){
			C_mat[r+MHS][r1].real( (-1.0)*l_coup*sin(PI/3) );
                        C_mat[r+MHS][r1].imag( (-1.0)*l_coup*cos(PI/3) );
                        C_mat[r1][r+MHS].real( (-1.0)*l_coup*sin(PI/3) );
                        C_mat[r1][r+MHS].imag( (1.0)*l_coup*cos(PI/3) );

                        C_mat[r][r1+MHS].real( (1.0)*l_coup*sin(PI/3)  );
                        C_mat[r][r1+MHS].imag( (-1.0)*l_coup*cos(PI/3) );
                        C_mat[r1+MHS][r].real( (1.0)*l_coup*sin(PI/3) );
                        C_mat[r1+MHS][r].imag( (1.0)*l_coup*cos(PI/3) );

			C_mat[r+MHS][r3].real( 0.0 );
                        C_mat[r+MHS][r3].imag( (1.0)*l_coup );
                        C_mat[r][r+MHS].real( 0.0 );
                        C_mat[r][r+MHS].imag( (-1.0)*l_coup );

                        C_mat[r][r3+MHS].real( 0.0 );
                        C_mat[r][r3+MHS].imag( (1.0)*l_coup );
                        C_mat[r3+MHS][r].real( 0.0 );
                        C_mat[r3+MHS][r].imag( (-1.0)*l_coup );
                }
                }
                }

		if(ap==1){
                if(ix!=Lx-1){
                if(iy==Ly-1){
                        C_mat[r1+MHS][r].real( (-1.0)*l_coup*sin(PI/3) );
                        C_mat[r1+MHS][r].imag( (-1.0)*l_coup*cos(PI/3) );
                        C_mat[r][r1+MHS].real( (-1.0)*l_coup*sin(PI/3) );
                        C_mat[r][r1+MHS].imag( (1.0)*l_coup*cos(PI/3) );

                        C_mat[r1][r+MHS].real( (1.0)*l_coup*sin(PI/3)  );
                        C_mat[r1][r+MHS].imag( (-1.0)*l_coup*cos(PI/3) );
                        C_mat[r+MHS][r1].real( (1.0)*l_coup*sin(PI/3) );
                        C_mat[r+MHS][r1].imag( (1.0)*l_coup*cos(PI/3) );

			C_mat[r2+MHS][r].real( (1.0)*l_coup*sin(PI/3) );
                        C_mat[r2+MHS][r].imag( (-1.0)*l_coup*cos(PI/3) );
                        C_mat[r][r2+MHS].real( (1.0)*l_coup*sin(PI/3) );
                        C_mat[r][r2+MHS].imag( (1.0)*l_coup*cos(PI/3) );

                        C_mat[r2][r+MHS].real( (-1.0)*l_coup*sin(PI/3)  );
                        C_mat[r2][r+MHS].imag( (-1.0)*l_coup*cos(PI/3) );
                        C_mat[r+MHS][r2].real( (-1.0)*l_coup*sin(PI/3) );
                        C_mat[r+MHS][r2].imag( (1.0)*l_coup*cos(PI/3) );
                }
                }
                }

                if(ap==2){
                if(ix!=Lx-1){
                if(iy==Ly-1){
                        C_mat[r+MHS][r1].real( (-1.0)*l_coup*sin(PI/3) );
                        C_mat[r+MHS][r1].imag( (-1.0)*l_coup*cos(PI/3) );
                        C_mat[r1][r+MHS].real( (-1.0)*l_coup*sin(PI/3) );
                        C_mat[r1][r+MHS].imag( (1.0)*l_coup*cos(PI/3) );

                        C_mat[r][r1+MHS].real( (1.0)*l_coup*sin(PI/3)  );
                        C_mat[r][r1+MHS].imag( (-1.0)*l_coup*cos(PI/3) );
                        C_mat[r1+MHS][r].real( (1.0)*l_coup*sin(PI/3) );
                        C_mat[r1+MHS][r].imag( (1.0)*l_coup*cos(PI/3) );

			C_mat[r+MHS][r2].real( (1.0)*l_coup*sin(PI/3) );
                        C_mat[r+MHS][r2].imag( (-1.0)*l_coup*cos(PI/3) );
                        C_mat[r2][r+MHS].real( (1.0)*l_coup*sin(PI/3) );
                        C_mat[r2][r+MHS].imag( (1.0)*l_coup*cos(PI/3) );

                        C_mat[r][r2+MHS].real( (-1.0)*l_coup*sin(PI/3)  );
                        C_mat[r][r2+MHS].imag( (-1.0)*l_coup*cos(PI/3) );
                        C_mat[r2+MHS][r].real( (-1.0)*l_coup*sin(PI/3) );
                        C_mat[r2+MHS][r].imag( (1.0)*l_coup*cos(PI/3) );
                }
                }
                }



                if(ap==1){
                if(ix==Lx-1){
                if(iy==Ly-1){
			C_mat[r1+MHS][r].real( (-1.0)*l_coup*sin(PI/3) );
                        C_mat[r1+MHS][r].imag( (-1.0)*l_coup*cos(PI/3) );
                        C_mat[r][r1+MHS].real( (-1.0)*l_coup*sin(PI/3) );
                        C_mat[r][r1+MHS].imag( (1.0)*l_coup*cos(PI/3) );

                        C_mat[r1][r+MHS].real( (1.0)*l_coup*sin(PI/3)  );
                        C_mat[r1][r+MHS].imag( (-1.0)*l_coup*cos(PI/3) );
                        C_mat[r+MHS][r1].real( (1.0)*l_coup*sin(PI/3) );
                        C_mat[r+MHS][r1].imag( (1.0)*l_coup*cos(PI/3) );
                }
                }
                }

		if(ap==2){
                if(ix==Lx-1){
                if(iy==Ly-1){
                        C_mat[r+MHS][r1].real( (-1.0)*l_coup*sin(PI/3) );
                        C_mat[r+MHS][r1].imag( (-1.0)*l_coup*cos(PI/3) );
                        C_mat[r1][r+MHS].real( (-1.0)*l_coup*sin(PI/3) );
                        C_mat[r1][r+MHS].imag( (1.0)*l_coup*cos(PI/3) );

                        C_mat[r][r1+MHS].real( (1.0)*l_coup*sin(PI/3)  );
                        C_mat[r][r1+MHS].imag( (-1.0)*l_coup*cos(PI/3) );
                        C_mat[r1+MHS][r].real( (1.0)*l_coup*sin(PI/3) );
                        C_mat[r1+MHS][r].imag( (1.0)*l_coup*cos(PI/3) );
		}
		}
		}
	}
	}
}


if(Periodic==true){
for(int sp=0;sp<2;sp++){
for(int ap=0;ap<3;ap++){
        for(int ix=0;ix<Lx;ix++){
        for(int iy=0;iy<Ly;iy++){

                r  = (ap   + 3*ix     + 3*Lx*iy );

                if(Half_Periodic==false){
                if(iy==Ly-1){
		r3 = (ap-1 + 3*ix );

                  if(ap==1){
			C_mat[r3+MHS][r].real( 0.0 );
                        C_mat[r3+MHS][r].imag( (1.0)*l_coup );
                        C_mat[r][r3+MHS].real( 0.0 );
                        C_mat[r][r3+MHS].imag( (-1.0)*l_coup );

                        C_mat[r3][r+MHS].real( 0.0 );
                        C_mat[r3][r+MHS].imag( (1.0)*l_coup );
                        C_mat[r+MHS][r3].real( 0.0 );
                        C_mat[r+MHS][r3].imag( (-1.0)*l_coup );
			
                  }

		  if(ap==2){
			C_mat[r+MHS][r3].real( 0.0 );
                        C_mat[r+MHS][r3].imag( (1.0)*l_coup );
                        C_mat[r][r+MHS].real( 0.0 );
                        C_mat[r][r+MHS].imag( (-1.0)*l_coup );

                        C_mat[r][r3+MHS].real( 0.0 );
                        C_mat[r][r3+MHS].imag( (1.0)*l_coup );
                        C_mat[r3+MHS][r].real( 0.0 );
                        C_mat[r3+MHS][r].imag( (-1.0)*l_coup );
                  }

                }
                }

                if(ix==Lx-1){
		r2 = (ap-1 + 3*Lx*iy );
                  if(ap==1){
			C_mat[r2+MHS][r].real( (1.0)*l_coup*sin(PI/3) );
                        C_mat[r2+MHS][r].imag( (-1.0)*l_coup*cos(PI/3) );
                        C_mat[r][r2+MHS].real( (1.0)*l_coup*sin(PI/3) );
                        C_mat[r][r2+MHS].imag( (1.0)*l_coup*cos(PI/3) );

                        C_mat[r2][r+MHS].real( (-1.0)*l_coup*sin(PI/3)  );
                        C_mat[r2][r+MHS].imag( (-1.0)*l_coup*cos(PI/3) );
                        C_mat[r+MHS][r2].real( (-1.0)*l_coup*sin(PI/3) );
                        C_mat[r+MHS][r2].imag( (1.0)*l_coup*cos(PI/3) );
                  }

		  if(ap==2){
                        C_mat[r+MHS][r2].real( (1.0)*l_coup*sin(PI/3) );
                        C_mat[r+MHS][r2].imag( (-1.0)*l_coup*cos(PI/3) );
                        C_mat[r2][r+MHS].real( (1.0)*l_coup*sin(PI/3) );
                        C_mat[r2][r+MHS].imag( (1.0)*l_coup*cos(PI/3) );

                        C_mat[r][r2+MHS].real( (-1.0)*l_coup*sin(PI/3)  );
                        C_mat[r][r2+MHS].imag( (-1.0)*l_coup*cos(PI/3) );
                        C_mat[r2+MHS][r].real( (-1.0)*l_coup*sin(PI/3) );
                        C_mat[r2+MHS][r].imag( (1.0)*l_coup*cos(PI/3) );
                  }
                }

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

sign_B1=B1_mag/abs(B1_mag);
sign_B2=B2_mag/abs(B2_mag);

connec_out<<sign_B1<<"	"<<sign_B2<<"	"<<sign_B1<<endl;

return 0;
}


