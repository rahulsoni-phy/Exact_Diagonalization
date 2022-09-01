#include "functions.h"
#include <iomanip>
#include <assert.h>
#define PI 3.14159265359

int main(int argc, char** argv){

int Lx, Ly, Mat_Size, Mat_Half_Size;
double hop, c_pot, l_coup, B_mag, B0_mag;
int k_point,k1_ind,k2_ind;

string connection_out;

//----------------------Inputs-------------------//
l_coup=atof(argv[1]);
B_mag=atof(argv[2]);
B0_mag=B_mag*(-2.0);

k1_ind=atoi(argv[3]);
k2_ind=atoi(argv[4]);

Lx=atoi(argv[5]);
Ly=atoi(argv[6]);

c_pot=0.6;
hop=1.0;
Mat_Size=6;
//Lx=24;
//Ly=24;

connection_out=argv[7];
ofstream connec_out(connection_out.c_str());

//-----------------------------------------------//
Mat_2_Complex_doub C_mat;
C_mat.resize(Mat_Size);

for(int i=0;i<Mat_Size;i++){
    C_mat[i].resize(Mat_Size);
    for(int j=0;j<Mat_Size;j++){
        C_mat[i][j].real(0.0);
        C_mat[i][j].imag(0.0);
    }
}

//B_mag
C_mat[0][0].real(-B0_mag);C_mat[0][0].imag(0.0);
C_mat[1][1].real(B0_mag);C_mat[1][1].imag(0.0);
C_mat[2][2].real(-B0_mag);C_mat[2][2].imag(0.0);
C_mat[3][3].real(B0_mag);C_mat[3][3].imag(0.0);
C_mat[4][4].real(-B_mag - c_pot);C_mat[4][4].imag(0.0);
C_mat[5][5].real(B_mag - c_pot);C_mat[5][5].imag(0.0);

double k1, k2;
	k1=(2.0*PI*k1_ind)/(1.0*Lx);
        k2=(2.0*PI*k2_ind)/(1.0*Ly);


//TBM
C_mat[0][4].real((-1.0)*( 1.0+cos(k1)+cos(k2) ) );C_mat[0][4].imag((1.0)*( sin(k1)+sin(k2) ) );
C_mat[1][5].real((-1.0)*( 1.0+cos(k1)+cos(k2) ) );C_mat[1][5].imag((1.0)*( sin(k1)+sin(k2) ) );

C_mat[2][4].real((-1.0)*( 1.0+cos(k1)+cos(k2) ) );C_mat[2][4].imag((-1.0)*( sin(k1)+sin(k2) ) );
C_mat[3][5].real((-1.0)*( 1.0+cos(k1)+cos(k2) ) );C_mat[3][5].imag((-1.0)*( sin(k1)+sin(k2) ) );

C_mat[4][0].real((-1.0)*( 1.0+cos(k1)+cos(k2) ) );C_mat[4][0].imag((-1.0)*( sin(k1)+sin(k2) ) );
C_mat[4][2].real((-1.0)*( 1.0+cos(k1)+cos(k2) ) );C_mat[4][2].imag((1.0)*( sin(k1)+sin(k2) ) );

C_mat[5][1].real((-1.0)*( 1.0+cos(k1)+cos(k2) ) );C_mat[5][1].imag((-1.0)*( sin(k1)+sin(k2) ) );
C_mat[5][3].real((-1.0)*( 1.0+cos(k1)+cos(k2) ) );C_mat[5][3].imag((1.0)*( sin(k1)+sin(k2) ) );

//R_SOC
C_mat[0][5].real( (-1.0*l_coup)*( sin(k1 + (2.0*PI/3.0) ) + sin(k2 + (4.0*PI/3.0) ) ) );
C_mat[0][5].imag( (-1.0*l_coup)*( 1.0 + cos(k1 + (2.0*PI/3.0) ) + cos(k2 + (4.0*PI/3.0) ) ) );

C_mat[1][4].real( (-1.0*l_coup)*( sin(k1 - (2.0*PI/3.0) ) + sin(k2 - (4.0*PI/3.0) ) ) );
C_mat[1][4].imag( (-1.0*l_coup)*( 1.0 + cos(k1 - (2.0*PI/3.0) ) + cos(k2 - (4.0*PI/3.0) ) ) );


C_mat[2][5].real( (1.0*l_coup)*( sin(k1 - (2.0*PI/3.0) ) + sin(k2 - (4.0*PI/3.0) ) ) );
C_mat[2][5].imag( (-1.0*l_coup)*( 1.0 + cos(k1 - (2.0*PI/3.0) ) + cos(k2 - (4.0*PI/3.0) ) ) );

C_mat[3][4].real( (1.0*l_coup)*( sin(k1 + (2.0*PI/3.0) ) + sin(k2 + (4.0*PI/3.0) ) ) );
C_mat[3][4].imag( (-1.0*l_coup)*( 1.0 + cos(k1 + (2.0*PI/3.0) ) + cos(k2 + (4.0*PI/3.0) ) ) );


C_mat[4][1].real( (-1.0*l_coup)*( sin(k1 - (2.0*PI/3.0) ) + sin(k2 - (4.0*PI/3.0) ) ) );
C_mat[4][1].imag( (1.0*l_coup)*( 1.0 + cos(k1 - (2.0*PI/3.0) ) + cos(k2 - (4.0*PI/3.0) ) ) );

C_mat[4][3].real( (1.0*l_coup)*( sin(k1 + (2.0*PI/3.0) ) + sin(k2 + (4.0*PI/3.0) ) ) );
C_mat[4][3].imag( (1.0*l_coup)*( 1.0 + cos(k1 + (2.0*PI/3.0) ) + cos(k2 + (4.0*PI/3.0) ) ) );


C_mat[5][0].real( (-1.0*l_coup)*( sin(k1 + (2.0*PI/3.0) ) + sin(k2 + (4.0*PI/3.0) ) ) );
C_mat[5][0].imag( (1.0*l_coup)*( 1.0 + cos(k1 + (2.0*PI/3.0) ) + cos(k2 + (4.0*PI/3.0) ) ) );

C_mat[5][2].real( (1.0*l_coup)*( sin(k1 - (2.0*PI/3.0) ) + sin(k2 - (4.0*PI/3.0) ) ) );
C_mat[5][2].imag( (1.0*l_coup)*( 1.0 + cos(k1 - (2.0*PI/3.0) ) + cos(k2 - (4.0*PI/3.0) ) ) );


for (int i=0;i<Mat_Size;i++){
    for (int j=0;j<Mat_Size;j++)
    {
    connec_out<<C_mat[i][j]<<" ";
    }
    connec_out<<endl;
}

//K-path
    Mat_1_intpair k_path;
    k_path.clear();
    pair_int temp_pair;

	int kx_i,ky_i;

    //Gamma to K
    kx_i=0;ky_i=0;
    while(kx_i< int((2.0*Lx)/3.0)){
        ky_i = 1*int(kx_i/2);
        temp_pair.first = kx_i;
        temp_pair.second = ky_i;
        k_path.push_back(temp_pair);
        kx_i++;
    }

    //K to M
    kx_i= int((2.0*Lx)/3.0);ky_i=-1*int(kx_i/2);
    while(ky_i < 0 ){
        kx_i= int(Lx/2) - int (ky_i/2);
        temp_pair.first = kx_i;
        temp_pair.second = -ky_i;
        k_path.push_back(temp_pair);
        ky_i++;
    }

    //M to Gamma
    kx_i= int (Lx/2); ky_i=0;
    while(kx_i >=0 ){
        temp_pair.first = kx_i;
        temp_pair.second = ky_i;
        k_path.push_back(temp_pair);
        kx_i--;
    }


return 0;
}
