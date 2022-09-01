#include "functions.h"
#include <iomanip>
#include <assert.h>
#define PI 3.14159265359

//string Evals_out="Eigenvalues.txt";
//ofstream Evals_file_out(Evals_out.c_str());


int main(int argc, char** argv){

string connection_in;
int k_point,k1_i,k2_i;

k_point=atoi(argv[1]);
k1_i=atoi(argv[2]);
k2_i=atoi(argv[3]);

//k1_i=atoi(argv[1]);
//k2_i=atoi(argv[2]);

connection_in=argv[4];
ifstream connec_in(connection_in.c_str());

//----------------------Inputs-------------------//
int H_size=6;
Mat_2_Complex_doub Hamiltonian;

//-----------------------------------------------//

Hamiltonian.resize(H_size);
for (int i=0;i<H_size;i++){
Hamiltonian[i].resize(H_size);
}

string tmp_str;
for (int i=0;i<H_size;i++){
    for (int j=0;j<H_size;j++){
    connec_in>>tmp_str;
    Hamiltonian[i][j]=reading_pair(tmp_str);
    }
}

//Getting Evals of the matrix
   Mat_2_Complex_doub Eigenvectors;
   Eigenvectors.resize(H_size);
   Mat_1_doub EVALS;
   EVALS.resize(H_size);


   for(int j=0;j<H_size;j++){
         Eigenvectors[j].resize(H_size);
        }
    int LDA=H_size;
    int info;
    /* Local arrays */
    double* eval = new double[H_size];

    lapack_complex_double* mat = (lapack_complex_double *) calloc(H_size*H_size, sizeof(lapack_complex_double));

    for(int i=0;i<H_size;i++){
        for(int j=0;j<=i;j++){

                mat[i*(H_size)+j] = Hamiltonian[i][j].real() + Hamiltonian[i][j].imag()*I; //On LAPTOP

        }
    }


    info=LAPACKE_zheev(LAPACK_ROW_MAJOR,  'V', 'L',  H_size, mat , LDA, eval);

    /* Check for convergence */
    if( info > 0 ) {
        cout<< "The LAPACKE_dsyev failed to diagonalize."<<endl;
    }


   for(int i=0;i<H_size;i++){
       EVALS[i]=eval[i];
        for(int j=0;j<H_size;j++){
        Eigenvectors[j][i].real(lapack_complex_double_real(mat[i*(H_size)+j]));
        Eigenvectors[j][i].imag(lapack_complex_double_imag(mat[i*(H_size)+j]));
	}
   }

	cout<<k_point<<"	"<<k1_i<<"	"<<k2_i<<"	"<<eval[0]<<"	"<<eval[1]<<"	"<<eval[2]<<"	"<<eval[3]<<"	"<<eval[4]<<"	"<<eval[5]<<endl;

//cout<<k1_i<<"      "<<eval[0]<<"   "<<eval[1]<<"   "<<eval[2]<<"   "<<eval[3]<<"   "<<eval[4]<<"   "<<eval[5]<<endl;


return 0;
}
