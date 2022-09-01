#include "functions.h"
#include <iomanip>
#include <assert.h>
#define PI 3.14159265359

using namespace std;

/*std::string to_string(int i){
   std::stringstream s;
   s << i;
   return s.str();
}
*/
double Fermi_(double en_, double mu_){
	double ffn, temp;
	temp=0.001;

	ffn=1.0/(1.0 + exp( (en_ - mu_)/temp ) );

return ffn;
}

//********************************************MAIN CODE BELOW*********************************************//
int main (int argc, char** argv){

int M, N, L_, Np_;
double LAMBDA_SOC,B_mag;

N=atoi(argv[1]);
M=atoi(argv[2]);
L_=3*M*N;

bool Doub_Occ=false;
bool CC_Specific_State=false;
bool RASHBA_CURRENTS=false;
bool Akxw=false;
bool SPECTRAL_FUNCT=false;
bool Current=false;
bool OP=false;
bool DOS=false;
bool LDOS_4x2=false;
bool LDOS=false;
bool OPT_COND=false;
bool Calculate_unk=false;

double ne_,ne_up, ne_down;
ne_=atof(argv[3]);
Np_= (int) (2*L_*ne_);

string connections_in;
connections_in=argv[4];
ifstream connec_in(connections_in.c_str());

string Evals_out="Eigenvalues.txt";
ofstream Evals_file_out(Evals_out.c_str());

int H_size;
Mat_2_Complex_doub Hamiltonian;
Mat_1_doub VAL;
VAL.resize(2);
double temp_double;

H_size=2*L_;
int MHS=H_size/2;

//cout<<"fine until here 1"<<endl;
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

for(int i=0;i<2;i++){
connec_in>>temp_double;
VAL[i]=temp_double;
}

LAMBDA_SOC=VAL[0];
B_mag=VAL[1];

//*********************************DIAGONALIZING HAMILTONIAN*******************************//

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

for(int i=0;i<H_size;i++){
Evals_file_out<<eval[i]<<endl;
}

double mu, epsilon, dmu_by_dN, N_el, mu_in, N_temp, Delta, dmu_by_dN_min;
dmu_by_dN=0.01*(EVALS[H_size-1]-EVALS[0])*(1.0/(1.0*H_size) );
dmu_by_dN_min = 0.0001*(EVALS[H_size-1]-EVALS[0])*(1.0/(1.0*H_size) );

epsilon=1e-2;
mu=EVALS[0];

N_temp=100000;

int iters =0;

while( abs(N_temp- 1.0*(Np_) )>epsilon){
N_temp=0;
for(int n=0;n<H_size;n++){
N_temp+=Fermi_(EVALS[n],mu);
}

//cout<<iters<<"   "<<setprecision(8)<<mu<<"	"<<N_temp<<endl;
mu = mu + dmu_by_dN*( (1.0*(Np_)) - N_temp);

iters++;
if((iters%1000)==0){
dmu_by_dN = (1000.0/(10.0*iters))*dmu_by_dN;
dmu_by_dN = max(dmu_by_dN, dmu_by_dN_min);
}

}

cout<<"Calculated mu = "<<mu<<endl;
cout<<"Final N = "<<N_temp<<endl;


//----------------------------------Edge Currents-----------------------------------------//
if(Current==true){

string Current_per_link1_out="Current_per_link.txt";
ofstream Current_link1_file_out(Current_per_link1_out.c_str());

string Current_per_link2_out="Current_per_link_test.txt";
ofstream Current_link2_file_out(Current_per_link2_out.c_str());

int m, mp;
m=0;mp=0;

for(int s=0;s<2;s++){
for(int ix=0;ix<N;ix++){
for(int iy=0;iy<M;iy++){
for(int a=0;a<3;a++){
        m=s*L_ + a + 3*ix + 3*N*iy;

for(int sp=0;sp<2;sp++){
for(int ixp=0;ixp<N;ixp++){
for(int iyp=0;iyp<M;iyp++){
for(int ap=0;ap<3;ap++){
        mp=sp*L_ + ap + 3*ixp + 3*N*iyp;

if(m<mp){
	if(Hamiltonian[m][mp].real()!=0 || Hamiltonian[m][mp].imag()!=0){
        complex<double>  Current1;
        Current1.real(0.0);
        Current1.imag(0.0);
        for(int alpha=0; alpha<H_size; alpha++){
Current1 = Current1 + ( (  (Hamiltonian[m][mp])* (conj(Eigenvectors[alpha][m]))*((Eigenvectors[alpha][mp])) )  - (  (Hamiltonian[mp][m])* ((conj(Eigenvectors[alpha][mp]))*((Eigenvectors[alpha][m])))  ) )*Fermi_(EVALS[alpha],mu);
        }

	Current_link1_file_out<<ix<<"	"<<iy<<"	"<<a<<"		"<<s<<"		"<<ixp<<"	"<<iyp<<"	"<<ap<<"	"<<sp<<"	"<<(-1.0)*Current1.imag()<<endl;
	Current_link2_file_out<<m<<"	"<<mp<<"	"<<(-1.0)*Current1.imag()<<endl;
        }
}

} } } }

} } } }

}
//-----------------------------------------------------------------------------------------//

if(Doub_Occ==true){
int MHS=H_size/2,m_;

Mat_1_Complex_doub DO_;
DO_.resize(MHS);
for(int i=0;i<MHS;i++){
DO_[i].real(0.0);DO_[i].imag(0.0);
}

for(int ix=0;ix<N;ix++){
for(int iy=0;iy<M;iy++){
for(int ap=0;ap<3;ap++){
        m_=ap + 3*ix + 3*N*iy;

        for(int n=0;n<H_size;n++){
	for(int p=0;p<H_size;p++){
	DO_[m_]=DO_[m_] + Fermi_(EVALS[n],mu)*Fermi_(EVALS[p],mu)*( ( conj(Eigenvectors[n][m_]) )*( Eigenvectors[n][m_] )*( conj(Eigenvectors[p][m_+MHS]) )*( Eigenvectors[p][m_+MHS] ) - ( conj(Eigenvectors[p][m_]) )*( Eigenvectors[n][m_] )*( conj(Eigenvectors[n][m_+MHS]) )*( Eigenvectors[p][m_+MHS] ) );
	}
	}
cout<<m_<<"	"<<DO_[m_].real()<<endl;
}
}
}

}


//-------------------------------------Order Parameter-------------------------------------//
if(OP==true){
int MHS=H_size/2;
int Nc=M*N;

Mat_1_Complex_doub CD_local_full,CD_total,CD_up,CD_dn,CD_up_puc,CD_dn_puc,CD_total_puc;  //Charge densities n_e(r,a)=\sum_{s,s'}c^{dag}_{r,a,s} delta_{s,s'}c_{r,a,s'}
Mat_1_Complex_doub SD_x,SD_y,SD_z,SD_z_puc,SD_x_puc,SD_y_puc;			//Spin densities n^{v}_s(r,a)=\sum_{s,s'}c^{dag}_{r,a,s} tau^{v}_{s,s'}c_{r,a,s'}

CD_local_full.resize(H_size);
CD_total.resize(MHS);CD_up.resize(MHS);CD_dn.resize(MHS);
SD_x.resize(MHS);SD_y.resize(MHS);SD_z.resize(MHS);

CD_total_puc.resize(Nc);CD_up_puc.resize(Nc);CD_dn_puc.resize(Nc);
SD_x_puc.resize(Nc);SD_y_puc.resize(Nc);SD_z_puc.resize(Nc);

for(int i=0;i<H_size;i++){
CD_local_full[i].real(0.0);CD_local_full[i].imag(0.0);
}

for(int i=0;i<MHS;i++){
CD_total[i].real(0.0);CD_total[i].imag(0.0);
CD_up[i].real(0.0);CD_up[i].imag(0.0);
CD_dn[i].real(0.0);CD_dn[i].imag(0.0);
SD_x[i].real(0.0);SD_x[i].imag(0.0);
SD_y[i].real(0.0);SD_y[i].imag(0.0);
SD_z[i].real(0.0);SD_z[i].imag(0.0);
}

for(int i=0;i<Nc;i++){
CD_total_puc[i].real(0.0);CD_total_puc[i].imag(0.0);
CD_up_puc[i].real(0.0);CD_up_puc[i].imag(0.0);
CD_dn_puc[i].real(0.0);CD_dn_puc[i].imag(0.0);
SD_x_puc[i].real(0.0);SD_x_puc[i].imag(0.0);
SD_y_puc[i].real(0.0);SD_y_puc[i].imag(0.0);
SD_z_puc[i].real(0.0);SD_z_puc[i].imag(0.0);
}

int m_,m_f;
complex<double> iota(0.0,1.0);

for(int ix=0;ix<N;ix++){
for(int iy=0;iy<M;iy++){
for(int ap=0;ap<3;ap++){
	m_=ap + 3*ix + 3*N*iy;

	for(int n=0;n<H_size;n++){
        CD_up[m_] = CD_up[m_] + (conj(Eigenvectors[n][m_]))*((Eigenvectors[n][m_]))*Fermi_(EVALS[n],mu);
	CD_dn[m_] = CD_dn[m_] + (conj(Eigenvectors[n][m_+MHS]))*((Eigenvectors[n][m_+MHS]))*Fermi_(EVALS[n],mu);
	SD_x[m_] = SD_x[m_] + 0.5*( (conj(Eigenvectors[n][m_]))*((Eigenvectors[n][m_+MHS])) + (conj(Eigenvectors[n][m_+MHS]))*((Eigenvectors[n][m_])) )*Fermi_(EVALS[n],mu);
	SD_y[m_] = SD_y[m_] + ( (-1.0)*iota )*0.5*( (conj(Eigenvectors[n][m_]))*((Eigenvectors[n][m_+MHS])) - (conj(Eigenvectors[n][m_+MHS]))*((Eigenvectors[n][m_])) )*Fermi_(EVALS[n],mu);	
	}
CD_total[m_] = CD_total[m_] + 1.0*(CD_up[m_] + CD_dn[m_]);
SD_z[m_] = SD_z[m_] + 0.5*(CD_up[m_] - CD_dn[m_]);

//---------------------
for(int s=0;s<2;s++){
	m_f=s*MHS + ap + 3*ix + 3*N*iy;
	for(int n=0;n<H_size;n++){
        CD_local_full[m_f] = CD_local_full[m_f] + (conj(Eigenvectors[n][m_f]))*((Eigenvectors[n][m_f]))*Fermi_(EVALS[n],mu);
	}
}
//---------------------

}
}
}

complex<double> Total_Magnetization;
Total_Magnetization.real(0.0);Total_Magnetization.imag(0.0);
for(int m=0;m<MHS;m++){
Total_Magnetization +=((1.0)/MHS)*SD_z[m];
}

cout<<LAMBDA_SOC<<"	"<<B_mag<<"	"<<Total_Magnetization<<endl;

complex<double> N_total_up,N_total_dn,Sz_total;
N_total_up.real(0.0);N_total_dn.real(0.0);Sz_total.real(0.0);
N_total_up.imag(0.0);N_total_dn.imag(0.0);Sz_total.imag(0.0);

string CDL_out="Charge_Density_local.txt";
ofstream CDL_file_out(CDL_out.c_str());

string CD_out="Charge_Density_up_dn_tot.txt";
ofstream CD_file_out(CD_out.c_str());

string SD_out="Spin_Density_x_y_z.txt";
ofstream SD_file_out(SD_out.c_str());

int m_f1,m_f2;

for(int m=0;m<MHS;m++){
	N_total_up += CD_up[m]; 
	N_total_dn += CD_dn[m];
	Sz_total += SD_z[m];

CD_file_out<<m<<"	"<<CD_up[m].real()<<"	"<<CD_dn[m].real()<<"	"<<CD_total[m].real()<<endl;
SD_file_out<<m<<"	"<<SD_x[m].real()<<"	"<<SD_y[m].real()<<"	"<<SD_z[m].real()<<"	"<<SD_x[m].real()*SD_x[m].real()<<"	"<<SD_y[m].real()*SD_y[m].real()<<"	"<<SD_z[m].real()*SD_z[m].real()<<"	"<<endl;

//-----------------------------------------//
	m_f1=m;
	m_f2=MHS+m;
CDL_file_out<<m<<"       "<<CD_local_full[m_f1].real()+CD_local_full[m_f2].real()<<endl;
//-----------------------------------------//
}

string CD_per_cell_out="Charge_Density_up_dn_tot_per_unit_cell.txt";
ofstream CD_per_cell_file_out(CD_per_cell_out.c_str());

string SD_per_cell_out="Spin_Density_x_y_z_per_unit_cell.txt";
ofstream SD_per_cell_file_out(SD_per_cell_out.c_str());


for(int iy=0;iy<M;iy++){
for(int ix=0;ix<N;ix++){
	m_f=ix + N*iy;
for(int ap=0;ap<3;ap++){
        m_=ap + 3*m_f;

	CD_up_puc[m_f] += CD_up[m_];  
	CD_dn_puc[m_f] += CD_dn[m_];
	CD_total_puc[m_f] += CD_total[m_];

	SD_x_puc[m_f] += SD_x[m_];
	SD_y_puc[m_f] += SD_y[m_];
	SD_z_puc[m_f] += SD_z[m_];
}

CD_per_cell_file_out<<m_f<<"	"<<CD_up_puc[m_f].real()<<"	"<<CD_dn_puc[m_f].real()<<"	"<<CD_total_puc[m_f].real()<<endl;
SD_per_cell_file_out<<m_f<<"    "<<SD_x_puc[m_f].real()<<"	"<<SD_y_puc[m_f].real()<<"	"<<SD_z_puc[m_f].real()<<endl;
}
}

///****************************************************************///
string CDL_specific="Charge_Density_for_specific_eigenvector.txt";
ofstream CDL_specific_file_out(CDL_specific.c_str());

int EV_spec=atoi(argv[5]);
//double k1_val=atof(argv[6]);
//double k2_val=atof(argv[6]);
complex<double> VAL;
complex<double> IOTA_COMPLEX(0.0,1.0);
complex<double> ZERO_COMPLEX(0.0,0.0);
VAL.real(0.0);VAL.imag(0.0);


for(int iy=0;iy<M;iy++){
for(int ix=0;ix<N;ix++){
for(int ap=0;ap<3;ap++){
        m_=ap + 3*ix + 3*N*iy;

//	VAL= (1.0/sqrt(N*1.0))*( exp(-IOTA_COMPLEX*k1_val*(PI*2.0*ix) ) )* (Eigenvectors[EV_spec][m_]+ Eigenvectors[EV_spec][m_+MHS]);
	VAL= (Eigenvectors[EV_spec][m_]+Eigenvectors[EV_spec][m_+MHS]);

// VAL= (2.0/(sqrt(N*1.0+1.0)*sqrt(M*1.0+1.0) ) )*( exp(IOTA_COMPLEX*k1_val*(2*PI*ix) )*exp(IOTA_COMPLEX*k2_val*(2*PI*iy) ) )* (Eigenvectors[EV_spec][m_]+ Eigenvectors[EV_spec][m_+MHS]);

CDL_specific_file_out<<m_<<"	"<<abs(VAL)*abs(VAL)<<endl;
}
}
}

double CD_per_cell;
int m;
string CD_EV_specific="Charge_Density_vs_r1_for_specific_eigenvector.txt";
ofstream file_CD_EV(CD_EV_specific.c_str());

for(int i=0;i<N;i++){
CD_per_cell=0.0;
for (int j=0;j<M;j++){
for (int a=0;a<3;a++){
m = a + 3*i + 3*N*j;

CD_per_cell +=(1.0/(M*1.0))*(conj(Eigenvectors[EV_spec][m]+Eigenvectors[EV_spec][m+MHS]) * (Eigenvectors[EV_spec][m]+Eigenvectors[EV_spec][m+MHS])).real();

}
}

file_CD_EV<<i<<"      "<<CD_per_cell<<endl;
}



}
//-----------------------------------------------------------------------------------------//

//--------------------------------------LDOS-----------------------------------------------//
if(LDOS==true){
double eta=0.005;
double one_by_PI_=1/PI;
double w_min=-0.6;
double w_max=0.6;
double dw=0.0025;

int w_size= (int) ( (w_max - w_min)/dw );
complex<double> iota(0.0,1.0);
complex<double> Zero_Complex(0.0,0.0);

Mat_1_Complex_doub Niw_edge,Niw_bulk;
Niw_edge.resize(w_size);
Niw_bulk.resize(w_size);

int site_edge,site_bulk;
site_edge=0;site_bulk=0;
double w;
w=0;

string base("_edge_and_bulk.txt");
string head("Local_Density_of_states_for_ap_");

for(int ap=0;ap<3;ap++){

string LDOS_out_(head + to_string(ap) + base);
ofstream LDOS_file_out(LDOS_out_.c_str());

for(int o=0;o<w_size;o++){
w = w_min + o*dw;

Niw_edge[o]=Zero_Complex;
Niw_bulk[o]=Zero_Complex;

for(int iy=0;iy<M;iy++){
for(int ix=0;ix<N;ix++){

	if(ix==0 || ix==N-1 || iy==M-1 || iy==0){
	site_edge = ap + 3*ix + 3*N*iy;
	
	for(int n=0; n<H_size; n++){
	Niw_edge[o] += (1.0/(2.0*(M+N)) )*( ( (conj(Eigenvectors[n][site_edge]))*((Eigenvectors[n][site_edge])) ) + ( (conj(Eigenvectors[n][MHS+site_edge]))*((Eigenvectors[n][MHS + site_edge])) ) )*(one_by_PI_)*( (eta) / ( (w-EVALS[n])*(w-EVALS[n]) + (eta*eta) ) );
	}

	}

	if(ix!=0 && ix!=N-1){
        if(iy!=0 && iy!=M-1){
        
	site_bulk = ap + 3*ix + 3*N*iy;
        for(int n=0; n<H_size; n++){
                Niw_bulk[o] += (1.0/(1.0*(M-2)*(N-2)) )*( (conj(Eigenvectors[n][site_bulk]))*((Eigenvectors[n][site_bulk])) + (conj(Eigenvectors[n][MHS+site_bulk]))*((Eigenvectors[n][MHS+site_bulk])) )*(one_by_PI_)*( (eta) / ( (w-EVALS[n])*(w-EVALS[n]) + (eta*eta) ) );
        }       
        
	}
        }


}
}

LDOS_file_out<<w<<"       "<<Niw_edge[o].real()<<"		"<<Niw_bulk[o].real()<<endl;
}

}
//---end---//

}
//-----------------------------------------------------------------------------------------//


//--------------------------------------DOS------------------------------------------------//
if(DOS==true){
double eta=0.01;
double value;
complex<double> value1,value2,value3;
double one_by_PI_=1/PI;
double w_min=-5.0;
double w_max=5.0;
double dw=0.002;

string DOS_out="Density_of_states.txt";
ofstream DOS_file_out(DOS_out.c_str());

double w=w_min;

int site_1,site_2,site_3;
while(w<=w_max){
value=0.0;
//value1.real(0.0);value2.real(0.0);value3.real(0.0);
//value1.imag(0.0);value2.imag(0.0);value3.imag(0.0);
//DOS_file_out<< w<<"     ";
DOS_file_out<< w<<"     ";

/*
for(int ix=0; ix<N; ix++){
for(int iy=0; iy<M; iy++){
site_1 = 0; //0 + 3*ix + 3*N*iy;
site_2 = 1; //1 + 3*ix + 3*N*iy;
site_3 = 2; //2  + 3*ix + 3*N*iy;

for(int n=0; n<H_size; n++){
value1=value1+(1.0/(2.0*L_))*(one_by_PI_)*((eta)/((w-EVALS[n])*(w-EVALS[n])+(eta*eta)))*( (conj(Eigenvectors[n][site_1]))*((Eigenvectors[n][site_1])) + (conj(Eigenvectors[n][MHS+site_1]))*((Eigenvectors[n][MHS+site_1])) );

value2=value2+(1.0/(2.0*L_))*(one_by_PI_)*((eta)/((w-EVALS[n])*(w-EVALS[n])+(eta*eta)))*( (conj(Eigenvectors[n][site_2]))*((Eigenvectors[n][site_2])) + (conj(Eigenvectors[n][MHS+site_2]))*((Eigenvectors[n][MHS+site_2])) );

value3=value3+(1.0/(2.0*L_))*(one_by_PI_)*((eta)/((w-EVALS[n])*(w-EVALS[n])+(eta*eta)))*( (conj(Eigenvectors[n][site_3]))*((Eigenvectors[n][site_3])) + (conj(Eigenvectors[n][MHS+site_3]))*((Eigenvectors[n][MHS+site_3])) );
}

}
}
*/

for(int i=0; i<H_size; i++){
value=value+(1.0/(2.0*L_))*(one_by_PI_)*((eta)/((w-EVALS[i])*(w-EVALS[i])+(eta*eta)));
}

DOS_file_out<<value<<endl;
//DOS_file_out<<value1.real()<<"     "<<value2.real()<<"  "<<value3.real()<<endl;
w=w+dw;

}
}
//----------------------------------------------------------------------------------------//
if(LDOS_4x2==true){
double eta=0.1;
complex<double> value;
double one_by_PI_=1/PI;
double w_min=-0.25;
double w_max=0.25;
double dw=0.05;

string DOS_out="Local_Density_of_states_4x2.txt";
ofstream DOS_file_out(DOS_out.c_str());

double w=w_min;

while(w<=w_max){
value.real(0.0);value.imag(0.0);
DOS_file_out<< w<<"     ";

for(int site=3; site<6; site++){
for(int n=0; n<H_size; n++){
value=value+(1.0/(2.0*3.0))*(one_by_PI_)*((eta)/((w-EVALS[n])*(w-EVALS[n])+(eta*eta)))*( (conj(Eigenvectors[n][site]))*((Eigenvectors[n][site])) + (conj(Eigenvectors[n][MHS+site]))*((Eigenvectors[n][MHS+site])) );
}
}


DOS_file_out<<value.real()<<"	"<<value.imag()<<endl;
w=w+dw;

}
}




//-------------------------------------A(k,w)---------------------------------------------//
if(SPECTRAL_FUNCT==true){
double w_min=-5;
double w_max=4.5;
double dw=0.08;
double one_by_L_=1.0/(2.0*L_);
double one_by_PI_=1.0/PI;
double eta=0.15;

int w_size= (int) ( (w_max - w_min)/dw );

complex<double> iota(0.0,1.0);

string file_out_SPECTRAL_fn = "SPECTRAL_FUNCTION.txt";
ofstream file_out_SPECTRAL(file_out_SPECTRAL_fn.c_str());

int a_size=3;
int s_size=2;

Mat_7_Complex_doub B_mat, A_mat;
B_mat.resize(a_size);
A_mat.resize(a_size);

for(int a=0; a<a_size; a++){
B_mat[a].resize(s_size);
for(int s=0; s<s_size; s++){
	B_mat[a][s].resize(a_size);
	for(int ap=0; ap<a_size; ap++){
		B_mat[a][s][ap].resize(s_size);
		for(int sp=0; sp<s_size; sp++){
			B_mat[a][s][ap][sp].resize(H_size);
			for(int j=0;j<H_size; j++){
				B_mat[a][s][ap][sp][j].resize(H_size);
				for(int jp=0;jp<H_size; jp++){
	                        	B_mat[a][s][ap][sp][j][jp].resize(w_size);
} } } } } }


for(int a=0; a<a_size; a++){
A_mat[a].resize(s_size);
for(int s=0; s<s_size; s++){
        A_mat[a][s].resize(a_size);
        for(int ap=0; ap<a_size; ap++){
                A_mat[a][s][ap].resize(s_size);
                for(int sp=0; sp<s_size; sp++){
                        A_mat[a][s][ap][sp].resize(N);
                        for(int k1=0;k1<N; k1++){
                                A_mat[a][s][ap][sp][k1].resize(M);
                                for(int k2=0;k2<M; k2++){
                                        A_mat[a][s][ap][sp][k1][k2].resize(w_size);
} } } } } }


for(int a=0; a<a_size; a++){
for(int s=0; s<s_size; s++){
for(int ap=0; ap<a_size; ap++){
for(int sp=0; sp<s_size; sp++){
for(int j=0;j<H_size; j++){
for(int jp=0;jp<H_size; jp++){
for(int o=0;o<w_size; o++){

B_mat[a][s][ap][sp][j][jp][o].real(0.0);
B_mat[a][s][ap][sp][j][jp][o].imag(0.0);

} } } } } } }


for(int a=0; a<a_size; a++){
for(int s=0; s<s_size; s++){
for(int ap=0; ap<a_size; ap++){
for(int sp=0; sp<s_size; sp++){
for(int k1=0;k1<N; k1++){
for(int k2=0;k2<M; k2++){
for(int o=0;o<w_size; o++){

A_mat[a][s][ap][sp][k1][k2][o].real(0.0);
A_mat[a][s][ap][sp][k1][k2][o].imag(0.0);

} } } } } } }


double w;
w=0.0;

int m,mp;
m=0;mp=0;

for(int s=0; s<s_size; s++){
for(int ix=0; ix<N; ix++){
for(int iy=0; iy<M; iy++){
for(int a=0; a<a_size; a++){
	m =  s*L_ + a + 3*ix + 3*iy*N;

for(int sp=0; sp<s_size; sp++){
for(int ixp=0; ixp<N; ixp++){
for(int iyp=0; iyp<M; iyp++){
for(int ap=0; ap<a_size; ap++){
        mp = sp*L_ + ap + 3*ixp + 3*iyp*N;

for(int o=0; o<w_size; o++){
	w=w_min+o*dw;

for(int n=0; n<H_size; n++){

B_mat[a][s][ap][sp][m][mp][o] += (one_by_PI_)*( ( conj(Eigenvectors[n][mp]) ) * (Eigenvectors[n][m]) )*( (eta)/((w-EVALS[n])*(w-EVALS[n])+(eta*eta)) );

} 

} 

} } } } 

} } } }

string base("_file.txt");
string head("Akw_a");

string ctr1("_s");
string ctr2("_ap");
string ctr3("_sp");

double k1,k2;

k1=0.0;
k2=0.0;

for(int s=0; s<s_size; s++){
for(int sp=s; sp<s_size; sp++){
for(int a=0; a<a_size; a++){
for(int ap=a; ap<a_size; ap++){

string file_new_(head + to_string(a) + ctr1 + to_string(s) + ctr2 + to_string(ap) + ctr3 + to_string(sp) + base);
ofstream file_new_in(file_new_.c_str());

for(int k1_ind=0; k1_ind<N; k1_ind++){
        k1 = (2.0 * PI * k1_ind) /(N*1.0) ;

for(int k2_ind=0; k2_ind<M; k2_ind++){
        k2 = (2.0 * PI * k2_ind) /(M*1.0) ;

for(int o=0; o<w_size; o++){
        w=w_min+o*dw;

for(int ix=0; ix<N; ix++){
for(int iy=0; iy<M; iy++){
        m =  s*L_ + a + 3*ix + 3*iy*N ;

for(int ixp=0; ixp<N; ixp++){
for(int iyp=0; iyp<M; iyp++){
        mp = sp*L_ + ap + 3*ixp + 3*iyp*N ;


A_mat[a][s][ap][sp][k1_ind][k2_ind][o] += (one_by_L_)*( exp( iota*k1*((ix-ixp)*1.0) )* exp( iota*( 1.0*k2*((iy-iyp)*1.0)) ) )*B_mat[a][s][ap][sp][m][mp][o]; 

} }

} }

file_new_in<<k1<<"      "<<k2<<"        "<<w<<"         "<<A_mat[a][s][ap][sp][k1_ind][k2_ind][o].real()<<endl;

}

} }

} } } }


complex<double> Spectral_Fn;
Spectral_Fn.real(0.0);
Spectral_Fn.imag(0.0);


string file_Spectral_p("Spectral_Function_path.txt");
ofstream file_Spectral_path(file_Spectral_p.c_str());

string file_Spectral_("Spectral_Function_final.txt");
ofstream file_Spectral_final(file_Spectral_.c_str());

for(int k1_ind=0; k1_ind<N; k1_ind++){
        k1 = (2.0 * PI * k1_ind) /(N*1.0) ;

for(int k2_ind=0; k2_ind<M; k2_ind++){
        k2 = (2.0 * PI * k2_ind) /(M*1.0) ;

for(int o=0; o<w_size; o++){
        w=w_min+o*dw;

for(int s=0; s<s_size; s++){
for(int a=0; a<a_size; a++){

Spectral_Fn += A_mat[a][s][a][s][k1_ind][k2_ind][o];

} }

file_Spectral_final<<k1<<"	"<<k2<<"	"<<w<<"		"<<Spectral_Fn.real()<<endl;

}

} }

cout << "Here 1"<<endl;

Mat_1_int K1_path, K2_path;
int k_index=0;
int n1, n2;
int dummy_n2;
assert(N%3 == 0);
assert(M%3 == 0);



//Gamma to K
cout <<"gamma to K"<<endl;
n1=0;n2=0;
while(n1< int((2.0*N)/3.0)){
n2 = -1*int(n1/2);

for(int o=0; o<w_size; o++){
        w=w_min+o*dw;
Spectral_Fn.real(0.0);
Spectral_Fn.imag(0.0);

for(int s=0; s<s_size; s++){
for(int a=0; a<a_size; a++){


dummy_n2=-n2;
Spectral_Fn += A_mat[a][s][a][s][n1][dummy_n2][o];


} }

file_Spectral_path<<k_index<<"   "<<n1<<"      "<<n2<<"        "<<w<<"         "<<Spectral_Fn.real()<<endl;

}

file_Spectral_path<<endl;
n1++;
k_index++;
cout << k_index<<endl;
}

//K to M
cout <<"K to m"<<endl;
n1= int((2.0*N)/3.0);n2=-1*int(n1/2);

while(n2 < 0 ){

n1= int(N/2) - int (n2/2);

for(int o=0; o<w_size; o++){
        w=w_min+o*dw;


Spectral_Fn.real(0.0);
Spectral_Fn.imag(0.0);

for(int s=0; s<s_size; s++){
for(int a=0; a<a_size; a++){

dummy_n2=-n2;
Spectral_Fn += A_mat[a][s][a][s][n1][dummy_n2][o];

} }

file_Spectral_path<<k_index<<"   "<<n1<<"      "<<n2<<"        "<<w<<"         "<<Spectral_Fn.real()<<endl;

}

file_Spectral_path<<endl;
n2++;
k_index++;
cout << k_index<<endl;
}


//M to Gamma
cout <<"m to gamma"<<endl;
n1= int (N/2); n2=0;
while(n1 >=0 ){

for(int o=0; o<w_size; o++){
        w=w_min+o*dw;
Spectral_Fn.real(0.0);
Spectral_Fn.imag(0.0);

for(int s=0; s<s_size; s++){
for(int a=0; a<a_size; a++){


Spectral_Fn += A_mat[a][s][a][s][n1][n2][o];

} }


file_Spectral_path<<k_index<<"   "<<n1<<"      "<<n2<<"        "<<w<<"         "<<Spectral_Fn.real()<<endl;
}

file_Spectral_path<<endl;
n1--;
k_index++;
cout << k_index<<endl;

}

}
//-----------------------------------------------------------------------------------------//

//---------------------------------------Akxw----------------------------------------------//
if(Akxw==true){
double w_min=-0.2;
double w_max=0.2;
double dw=0.001;
double one_by_L_=1.0/(2.0*L_);
double one_by_PI_=1.0/PI;
double eta=0.002;

int w_size= (int) ( (w_max - w_min)/dw );
complex<double> iota(0.0,1.0);

int a_size=3;
int s_size=2;

Mat_3_Complex_doub B_mat;
B_mat.resize(N);

for(int i=0;i<N;i++){
B_mat[i].resize(N);
	for(int j=0;j<N;j++){
	B_mat[i][j].resize(w_size);
	}
}

cout << "r1"<<endl;

for(int i=0;i<N;i++){
	for(int j=0;j<N;j++){
		for(int o=0;o<w_size;o++){
		B_mat[i][j][o].real(0.0);
		B_mat[i][j][o].imag(0.0);
		}
	}
}

double w;
w=0.0;

int m,mp;
m=0;mp=0;

for(int ix=0; ix<N; ix++){
for(int ixp=0; ixp<N; ixp++){

for(int o=0; o<w_size; o++){
	w=w_min+o*dw;

for(int iy=0; iy<M; iy++){
for(int s=0; s<s_size; s++){
for(int a=0; a<a_size; a++){
	m =  s*L_ + a + 3*ix + 3*iy*N;
	mp = s*L_ + a + 3*ixp + 3*iy*N;

for(int n=0; n<H_size; n++){

B_mat[ix][ixp][o] += (one_by_PI_)*( ( conj(Eigenvectors[n][mp]) ) * (Eigenvectors[n][m]) )*( (eta)/((w-EVALS[n])*(w-EVALS[n])+(eta*eta)) );

}

}}

} 

} 

}}

double k1;
k1=0.0;

complex<double> SP_fn;

string file_out_Akxw_ = "Akxw.txt";
ofstream file_out_Akxw(file_out_Akxw_.c_str());

for(int k1_ind=0; k1_ind<N; k1_ind++){
k1 = (2.0 * PI * k1_ind) /(N*1.0);

	for(int o=0; o<w_size; o++){
	SP_fn.real(0.0);SP_fn.imag(0.0);
	w=w_min+o*dw;

		for(int ix=0; ix<N; ix++){
		for(int ixp=0; ixp<N; ixp++){
		
		SP_fn += exp(iota*(k1*(ix - ixp)*1.0))*B_mat[ix][ixp][o];

		}	
		}
	file_out_Akxw<<k1_ind<<"	"<<w<<"		"<<SP_fn.real()<<"	"<<SP_fn.imag()<<endl;
	}

file_out_Akxw<<endl;
}

}
//-----------------------------------------------------------------------------------------//

//---------------------------------------Akxw----------------------------------------------//
if(Akxw==true){
double w_min=-0.2;
double w_max=0.2;
double dw=0.001;
double one_by_L_=1.0/(2.0*L_);
double one_by_PI_=1.0/PI;
double eta=0.002;

int w_size= (int) ( (w_max - w_min)/dw );
complex<double> iota(0.0,1.0);

int a_size=3;
int s_size=2;

Mat_3_Complex_doub B_mat;
B_mat.resize(N);

for(int i=0;i<N;i++){
B_mat[i].resize(N);
        for(int j=0;j<N;j++){
        B_mat[i][j].resize(w_size);
        }
}

cout << "r1"<<endl;

for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
                for(int o=0;o<w_size;o++){
                B_mat[i][j][o].real(0.0);
                B_mat[i][j][o].imag(0.0);
                }
        }
}

double w;
w=0.0;

int m,mp;
m=0;mp=0;

for(int ix=0; ix<N; ix++){
for(int ixp=0; ixp<N; ixp++){

for(int o=0; o<w_size; o++){
        w=w_min+o*dw;

for(int iy=0; iy<M; iy++){
for(int s=0; s<s_size; s++){
for(int a=0; a<a_size; a++){
        m =  s*L_ + a + 3*ix + 3*iy*N;
        mp = s*L_ + a + 3*ixp + 3*iy*N;

for(int n=0; n<H_size; n++){

B_mat[ix][ixp][o] += (one_by_PI_)*( ( conj(Eigenvectors[n][mp]) ) * (Eigenvectors[n][m]) )*( (eta)/((w-EVALS[n])*(w-EVALS[n])+(eta*eta)) );

}

}}

}

}

}}

double k1,k2;
k1=0.0;k2=0.0;

complex<double> SP_fn;

string file_out_Akxw_ = "Akxw_OxO.txt";
ofstream file_out_Akxw(file_out_Akxw_.c_str());

for(int k1_ind=1; k1_ind<N; k1_ind++){
k1 = ( PI * k1_ind) /(N*1.0+1.0);

        for(int o=0; o<w_size; o++){
        SP_fn.real(0.0);SP_fn.imag(0.0);
        w=w_min+o*dw;

	for(int k2_ind=1; k2_ind<N; k2_ind++){
		k2 = ( PI * k2_ind) /(M*1.0+1.0);
	
                for(int ix=0; ix<N; ix++){
                for(int ixp=0; ixp<N; ixp++){

                SP_fn += sin((k1*ix*1.0)*sin((k1*ixp)*1.0))*B_mat[ix][ixp][o];

                }
                }
	}

        file_out_Akxw<<k1_ind<<"        "<<w<<"         "<<SP_fn.real()<<"      "<<SP_fn.imag()<<endl;
        }

file_out_Akxw<<endl;
}

}
//-----------------------------------------------------------------------------------------//




//------------------------------------Proper Currents--------------------------------------//
if(RASHBA_CURRENTS==true){

complex<double> Current_temp1,Current_temp2;
Current_temp1.real(0.0);Current_temp1.imag(0.0);
Current_temp2.real(0.0);Current_temp2.imag(0.0);

complex<double> temp_val1,temp_val2;
temp_val1.real(0.0);temp_val1.imag(0.0);
temp_val2.real(0.0);temp_val2.imag(0.0);

int a_size,s_size;
a_size=3;s_size=2;

Mat_4_Complex_doub CURRENT;
Mat_2_Complex_doub CHARGE_CURRENT_UP_UP,CHARGE_CURRENT_DN_DN,CHARGE_CURRENT_UP_DN,CHARGE_CURRENT_DN_UP;

Mat_2_Complex_doub SPIN_CURRENT_UP_UP,SPIN_CURRENT_DN_DN,SPIN_CURRENT_UP_DN,SPIN_CURRENT_DN_UP;

int Lattice_size=3*M*N;

CURRENT.resize(Lattice_size);
for (int i=0;i<Lattice_size;i++){
CURRENT[i].resize(Lattice_size);
	for(int j=0;j<Lattice_size;j++){
	CURRENT[i][j].resize(s_size);
		for (int s=0;s<s_size;s++){
		CURRENT[i][j][s].resize(s_size);
		}
	}
}

int m,mp,temp_m,temp_mp;
complex<double> iota(0.0,1.0);
complex<double> Zero_Complex(0.0,0.0);

string CDC_Matrix="CDC_UU_DD_UD_DU.txt";
ofstream file_CDC_Matrix(CDC_Matrix.c_str());

for(int i1=0;i1<N;i1++){
for(int i2=0;i2<M;i2++){
for(int a=0;a<a_size;a++){
	m = a + 3*i1 + 3*N*i2;
	for(int i1p=0;i1p<N;i1p++){
	for(int i2p=0;i2p<M;i2p++){
	for(int ap=0;ap<a_size;ap++){
		mp = ap + 3*i1p + 3*N*i2p;
		for(int s=0;s<s_size;s++){
		temp_m = s*L_ + m;
		for(int sp=0;sp<s_size;sp++){
		temp_mp = sp*L_ + mp;
		
                        CURRENT[m][mp][s][sp]=Zero_Complex;
			Current_temp1=Zero_Complex;

			if(Hamiltonian[temp_m][temp_mp].real()!=0 || Hamiltonian[temp_m][temp_mp].imag()!=0){
				for(int n=0;n<H_size;n++){
				Current_temp1 += ( (conj(Eigenvectors[n][temp_m]))*((Eigenvectors[n][temp_mp])) )*Fermi_(EVALS[n],mu);
                                }
                        CURRENT[m][mp][s][sp] = Current_temp1;
			}

		}
		}
	}
	}
	}
}
}
}

for(int q=0;q<Lattice_size;q++){
	for(int qp=0;qp<Lattice_size;qp++){
	if(abs(CURRENT[q][qp][0][0] + CURRENT[q][qp][1][1] + CURRENT[q][qp][0][1] + CURRENT[q][qp][1][0])!=0){
	file_CDC_Matrix<<q<<"      "<<qp<<"        "<<CURRENT[q][qp][0][0]<<"      "<<CURRENT[q][qp][1][1]<<"        "<<CURRENT[q][qp][0][1]<<"      "<<CURRENT[q][qp][1][0]<<endl;
	}
	}
}

CHARGE_CURRENT_UP_UP.resize(Lattice_size);
CHARGE_CURRENT_DN_DN.resize(Lattice_size);
CHARGE_CURRENT_UP_DN.resize(Lattice_size);
CHARGE_CURRENT_DN_UP.resize(Lattice_size);
for (int i=0;i<Lattice_size;i++){
CHARGE_CURRENT_UP_UP[i].resize(Lattice_size);
CHARGE_CURRENT_DN_DN[i].resize(Lattice_size);
CHARGE_CURRENT_UP_DN[i].resize(Lattice_size);
CHARGE_CURRENT_DN_UP[i].resize(Lattice_size);
}

SPIN_CURRENT_UP_UP.resize(Lattice_size);
SPIN_CURRENT_DN_DN.resize(Lattice_size);
SPIN_CURRENT_UP_DN.resize(Lattice_size);
SPIN_CURRENT_DN_UP.resize(Lattice_size);
for (int i=0;i<Lattice_size;i++){
SPIN_CURRENT_UP_UP[i].resize(Lattice_size);
SPIN_CURRENT_DN_DN[i].resize(Lattice_size);
SPIN_CURRENT_UP_DN[i].resize(Lattice_size);
SPIN_CURRENT_DN_UP[i].resize(Lattice_size);
}


for (int i=0;i<N;i++){
for(int j=0;j<M;j++){
for (int a=0;a<a_size;a++){
	m = a + 3*i + 3*N*j;
	for (int ip=0;ip<N;ip++){
	for(int jp=0;jp<M;jp++){
	for (int ap=0;ap<a_size;ap++){
		mp = ap + 3*ip + 3*N*jp;
	CHARGE_CURRENT_UP_UP[m][mp]=Zero_Complex;
	CHARGE_CURRENT_DN_DN[m][mp]=Zero_Complex;
	CHARGE_CURRENT_UP_DN[m][mp]=Zero_Complex;
	CHARGE_CURRENT_DN_UP[m][mp]=Zero_Complex;

	SPIN_CURRENT_UP_UP[m][mp]=Zero_Complex;
        SPIN_CURRENT_DN_DN[m][mp]=Zero_Complex;
        SPIN_CURRENT_UP_DN[m][mp]=Zero_Complex;
        SPIN_CURRENT_DN_UP[m][mp]=Zero_Complex;
	}
	}
	}
}
}
}


for (int i=0;i<N;i++){
for(int j=0;j<M;j++){
for (int a=0;a<a_size;a++){
m = a + 3*i + 3*N*j;
        for (int ip=0;ip<N;ip++){
        for(int jp=0;jp<M;jp++){
        for (int ap=0;ap<a_size;ap++){
	mp = ap + 3*ip + 3*N*jp;

	CHARGE_CURRENT_UP_UP[m][mp] = -iota*(CURRENT[m][mp][0][0] - CURRENT[mp][m][0][0]);
	CHARGE_CURRENT_DN_DN[m][mp] = -iota*(CURRENT[m][mp][1][1] - CURRENT[mp][m][1][1]);

//	CHARGE_CURRENT_UP_UP[m][mp] = -iota*(CURRENT[m][mp][0][0] - CURRENT[mp][m][0][0]);
//	CHARGE_CURRENT_DN_DN[m][mp] = -iota*(CURRENT[m][mp][1][1] - CURRENT[mp][m][1][1]);	
//	CHARGE_CURRENT_UP_DN[m][mp] = 0.15*(CURRENT[m][mp][0][1] +  CURRENT[mp][m][1][0] + CURRENT[m][mp][1][0] + CURRENT[mp][m][0][1]  );
//	CHARGE_CURRENT_DN_UP[m][mp] = 0.15*(CURRENT[m][mp][0][1] +  CURRENT[mp][m][1][0] - CURRENT[m][mp][1][0] - CURRENT[mp][m][0][1]  );
	

	if(ip==i && jp==j){
		if(a==1){
			if(ap==0){
			CHARGE_CURRENT_UP_DN[mp][m] = (1.0)*LAMBDA_SOC*(CURRENT[mp][m][0][1] + CURRENT[m][mp][1][0]);
			CHARGE_CURRENT_DN_UP[mp][m] = (1.0)*LAMBDA_SOC*(CURRENT[mp][m][1][0] + CURRENT[m][mp][0][1]);

			CHARGE_CURRENT_UP_DN[m][mp] = (-1.0)*LAMBDA_SOC*(CURRENT[m][mp][0][1] + CURRENT[mp][m][1][0]);
			CHARGE_CURRENT_DN_UP[m][mp] = (-1.0)*LAMBDA_SOC*(CURRENT[m][mp][1][0] + CURRENT[mp][m][0][1]);
			}
	
			if(ap==2){
			CHARGE_CURRENT_UP_DN[mp][m] = (1.0)*LAMBDA_SOC*(CURRENT[mp][m][0][1] + CURRENT[m][mp][1][0]);
			CHARGE_CURRENT_DN_UP[mp][m] = (1.0)*LAMBDA_SOC*(CURRENT[mp][m][1][0] + CURRENT[m][mp][0][1]);

			CHARGE_CURRENT_UP_DN[m][mp] = (-1.0)*LAMBDA_SOC*(CURRENT[m][mp][0][1] + CURRENT[mp][m][1][0]);
			CHARGE_CURRENT_DN_UP[m][mp] = (-1.0)*LAMBDA_SOC*(CURRENT[m][mp][1][0] + CURRENT[mp][m][0][1]);
			}
		}
	}

	if(jp==j-1 && j!=0){
		if(a==1){
			if(ap==2){
			CHARGE_CURRENT_UP_DN[mp][m] = (1.0)*LAMBDA_SOC*(exp((2.0)*iota*PI/3.0)*CURRENT[mp][m][0][1] + exp((-2.0)*iota*PI/3.0)*CURRENT[m][mp][1][0]);
			CHARGE_CURRENT_DN_UP[mp][m] = (1.0)*LAMBDA_SOC*(exp((-2.0)*iota*PI/3.0)*CURRENT[mp][m][1][0] + exp((2.0)*iota*PI/3.0)*CURRENT[m][mp][0][1]);

			CHARGE_CURRENT_UP_DN[m][mp] = (-1.0)*LAMBDA_SOC*(exp((2.0)*iota*PI/3.0)*CURRENT[m][mp][0][1] + exp((-2.0)*iota*PI/3.0)*CURRENT[mp][m][1][0]);
			CHARGE_CURRENT_DN_UP[m][mp] = (-1.0)*LAMBDA_SOC*(exp((-2.0)*iota*PI/3.0)*CURRENT[m][mp][1][0] + exp((2.0)*iota*PI/3.0)*CURRENT[mp][m][0][1]);
			}
		}
	}

	if(jp==j+1 && j!=M-1){
		if(a==1){
			if(ap==0){
			CHARGE_CURRENT_UP_DN[mp][m] = (1.0)*LAMBDA_SOC*(exp((2.0)*iota*PI/3.0)*CURRENT[mp][m][0][1] + exp((-2.0)*iota*PI/3.0)*CURRENT[m][mp][1][0]);
			CHARGE_CURRENT_DN_UP[mp][m] = (1.0)*LAMBDA_SOC*(exp((-2.0)*iota*PI/3.0)*CURRENT[mp][m][1][0] + exp((2.0)*iota*PI/3.0)*CURRENT[m][mp][0][1]);

			CHARGE_CURRENT_UP_DN[m][mp] = (-1.0)*LAMBDA_SOC*(exp((2.0)*iota*PI/3.0)*CURRENT[m][mp][0][1] + exp((-2.0)*iota*PI/3.0)*CURRENT[mp][m][1][0]);
			CHARGE_CURRENT_DN_UP[m][mp] = (-1.0)*LAMBDA_SOC*(exp((-2.0)*iota*PI/3.0)*CURRENT[m][mp][1][0] + exp((2.0)*iota*PI/3.0)*CURRENT[mp][m][0][1]);
			}
		}
	}

	if(ip==i+1 && i!=N-1){
		if(a==1){
			if(ap==0){
			CHARGE_CURRENT_UP_DN[mp][m] = (1.0)*LAMBDA_SOC*(exp((-2.0)*iota*PI/3.0)*CURRENT[mp][m][0][1] + exp((2.0)*iota*PI/3.0)*CURRENT[m][mp][1][0]);
			CHARGE_CURRENT_DN_UP[mp][m] = (1.0)*LAMBDA_SOC*(exp((2.0)*iota*PI/3.0)*CURRENT[mp][m][1][0] + exp((-2.0)*iota*PI/3.0)*CURRENT[m][mp][0][1]);

			CHARGE_CURRENT_UP_DN[m][mp] = (-1.0)*LAMBDA_SOC*(exp((-2.0)*iota*PI/3.0)*CURRENT[m][mp][0][1] + exp((2.0)*iota*PI/3.0)*CURRENT[mp][m][1][0]);
			CHARGE_CURRENT_DN_UP[m][mp] = (-1.0)*LAMBDA_SOC*(exp((2.0)*iota*PI/3.0)*CURRENT[m][mp][1][0] + exp((-2.0)*iota*PI/3.0)*CURRENT[mp][m][0][1]);
			}
		}
	}

	if(ip==i-1 && i!=0){
		if(a==1){
			if(ap==2){
			CHARGE_CURRENT_UP_DN[mp][m] = (1.0)*LAMBDA_SOC*(exp((-2.0)*iota*PI/3.0)*CURRENT[mp][m][0][1] + exp((2.0)*iota*PI/3.0)*CURRENT[m][mp][1][0]);
			CHARGE_CURRENT_DN_UP[mp][m] = (1.0)*LAMBDA_SOC*(exp((2.0)*iota*PI/3.0)*CURRENT[mp][m][1][0] + exp((-2.0)*iota*PI/3.0)*CURRENT[m][mp][0][1]);

			CHARGE_CURRENT_UP_DN[m][mp] = (-1.0)*LAMBDA_SOC*(exp((-2.0)*iota*PI/3.0)*CURRENT[m][mp][0][1] + exp((2.0)*iota*PI/3.0)*CURRENT[mp][m][1][0]);
			CHARGE_CURRENT_DN_UP[m][mp] = (-1.0)*LAMBDA_SOC*(exp((2.0)*iota*PI/3.0)*CURRENT[m][mp][1][0] + exp((-2.0)*iota*PI/3.0)*CURRENT[mp][m][0][1]);
			}
		}
	}

//-------------------Periodic-y currents----------------------------//
/*	if(jp==0 && j==M-1){
		if(a==1){
			if(ap==0){
			CHARGE_CURRENT_UP_DN[mp][m] = (1.0)*LAMBDA_SOC*(exp((2.0)*iota*PI/3.0)*CURRENT[mp][m][0][1] + exp((-2.0)*iota*PI/3.0)*CURRENT[m][mp][1][0]);
                        CHARGE_CURRENT_DN_UP[mp][m] = (1.0)*LAMBDA_SOC*(exp((-2.0)*iota*PI/3.0)*CURRENT[mp][m][1][0] + exp((2.0)*iota*PI/3.0)*CURRENT[m][mp][0][1]);

                        CHARGE_CURRENT_UP_DN[m][mp] = (-1.0)*LAMBDA_SOC*(exp((2.0)*iota*PI/3.0)*CURRENT[m][mp][0][1] + exp((-2.0)*iota*PI/3.0)*CURRENT[mp][m][1][0]);
                        CHARGE_CURRENT_DN_UP[m][mp] = (-1.0)*LAMBDA_SOC*(exp((-2.0)*iota*PI/3.0)*CURRENT[m][mp][1][0] + exp((2.0)*iota*PI/3.0)*CURRENT[mp][m][0][1]);
			}
		}
	}

	if(jp==M-1 && j==0){
                if(a==1){
                        if(ap==2){
			CHARGE_CURRENT_UP_DN[mp][m] = (1.0)*LAMBDA_SOC*(exp((2.0)*iota*PI/3.0)*CURRENT[mp][m][0][1] + exp((-2.0)*iota*PI/3.0)*CURRENT[m][mp][1][0]);
                        CHARGE_CURRENT_DN_UP[mp][m] = (1.0)*LAMBDA_SOC*(exp((-2.0)*iota*PI/3.0)*CURRENT[mp][m][1][0] + exp((2.0)*iota*PI/3.0)*CURRENT[m][mp][0][1]);

                        CHARGE_CURRENT_UP_DN[m][mp] = (-1.0)*LAMBDA_SOC*(exp((2.0)*iota*PI/3.0)*CURRENT[m][mp][0][1] + exp((-2.0)*iota*PI/3.0)*CURRENT[mp][m][1][0]);
                        CHARGE_CURRENT_DN_UP[m][mp] = (-1.0)*LAMBDA_SOC*(exp((-2.0)*iota*PI/3.0)*CURRENT[m][mp][1][0] + exp((2.0)*iota*PI/3.0)*CURRENT[mp][m][0][1]);
                        }
                }
	}
//-----------------------------------------------------------------//
*/	}
	}
	}
}
}
}

string CC_FULL="Total_DN_Charge_Current_at_each_link.txt";
ofstream file_CC_FULL(CC_FULL.c_str());

string CC_Individual="Individual_Charge_Current_at_each_link.txt";
ofstream file_CC_Individual(CC_Individual.c_str());

string SC_FULL="Average_Charge_Current_per_unit_cell_in_periodic_y.txt";
ofstream file_SC_FULL(SC_FULL.c_str());

file_CC_Individual<<"#  "<<"q"<<"      "<<"qp"<<"        "<<"UP->UP"<<"         "<<"DN->DN"<<"          "<<"DN->UP"<<"          "<<"UP->DN"<<endl;

for(int q=0;q<Lattice_size;q++){
	for(int qp=0;qp<Lattice_size;qp++){
	if(abs(CHARGE_CURRENT_UP_UP[q][qp])!=0 || abs(CHARGE_CURRENT_DN_DN[q][qp])!=0 || abs(CHARGE_CURRENT_UP_DN[q][qp])!=0 || abs(CHARGE_CURRENT_DN_UP[q][qp])!=0){
		if(q<qp){
//		file_CC_FULL<<q<<"      "<<qp<<"        "<<(CHARGE_CURRENT_UP_UP[qp][q] + CHARGE_CURRENT_DN_DN[qp][q] + CHARGE_CURRENT_UP_DN[qp][q] + CHARGE_CURRENT_DN_UP[qp][q]).real()<<endl;
		file_CC_FULL<<q<<"      "<<qp<<"        "<<(CHARGE_CURRENT_DN_UP[qp][q] + CHARGE_CURRENT_DN_DN[qp][q]).real()<<endl;
//		file_SC_FULL<<q<<"      "<<qp<<"        "<<( (CHARGE_CURRENT_UP_UP[qp][q] - CHARGE_CURRENT_DN_DN[qp][q] + CHARGE_CURRENT_DN_UP[qp][q] - CHARGE_CURRENT_UP_DN[qp][q]).real() )*0.5<<endl;
//		file_SC_FULL<<q<<"      "<<qp<<"        "<<( (CHARGE_CURRENT_UP_UP[qp][q] - CHARGE_CURRENT_DN_DN[qp][q] ).real() )<<endl;
		}

	file_CC_Individual<<q<<"	"<<qp<<"	"<<CHARGE_CURRENT_UP_UP[qp][q].real()<<"	"<<CHARGE_CURRENT_DN_DN[qp][q].real()<<"	"<<CHARGE_CURRENT_UP_DN[qp][q].real()<<"	"<<CHARGE_CURRENT_DN_UP[qp][q].real()<<endl;

	}
	}
}


complex<double> CC_per_cell;

for(int i=0;i<N;i++){
CC_per_cell=Zero_Complex;
for (int j=0;j<M;j++){
for (int a=0;a<a_size;a++){
m = a + 3*i + 3*N*j;
        for (int jp=0;jp<M;jp++){
        for (int ap=0;ap<a_size;ap++){
        mp = ap + 3*i + 3*N*jp;

	if(jp==j && ap==a+1){
	CC_per_cell+= (1.0/(M*1.0))*(CHARGE_CURRENT_UP_UP[mp][m] + CHARGE_CURRENT_DN_DN[mp][m] + CHARGE_CURRENT_UP_DN[mp][m] + CHARGE_CURRENT_DN_UP[mp][m]);
	}
	
	if(jp==j+1 && ap==a-1){
	CC_per_cell+= (1.0/(M*1.0))*(CHARGE_CURRENT_UP_UP[mp][m] + CHARGE_CURRENT_DN_DN[mp][m] + CHARGE_CURRENT_UP_DN[mp][m] + CHARGE_CURRENT_DN_UP[mp][m]);
	}

	if(jp==0 && j==M-1 && ap==a-1){
	CC_per_cell+= (1.0/(M*1.0))*(CHARGE_CURRENT_UP_UP[mp][m] + CHARGE_CURRENT_DN_DN[mp][m] + CHARGE_CURRENT_UP_DN[mp][m] + CHARGE_CURRENT_DN_UP[mp][m]);
	}
	
	}
	}
}
}

file_SC_FULL<<i<<"	"<<CC_per_cell.real()<<endl;
}



string CC_Conservation="Charge_Current_Conservation.txt";
ofstream CC_Conservation_out(CC_Conservation.c_str());

complex<double> Total_Charge_Current_at_site,Total_Spin_Current_at_site;

for(int j=0;j<M;j++){
	for(int i=0;i<N;i++){
		for(int a=0;a<a_size;a++){
		m = a + 3*i + 3*N*j;
		Total_Charge_Current_at_site.real(0.0);Total_Spin_Current_at_site.real(0.0);
		Total_Charge_Current_at_site.imag(0.0);Total_Spin_Current_at_site.imag(0.0);
		

		for(int ip=0;ip<N;ip++){
			for(int jp=0;jp<M;jp++){
				for(int ap=0;ap<a_size;ap++){
				mp = ap + 3*ip + 3*N*jp;
				if(abs(CHARGE_CURRENT_UP_UP[m][mp])!=0 || abs(CHARGE_CURRENT_DN_DN[m][mp])!=0 || abs(CHARGE_CURRENT_UP_DN[m][mp])!=0 || abs(CHARGE_CURRENT_DN_UP[m][mp])!=0){
				Total_Charge_Current_at_site += CHARGE_CURRENT_UP_UP[m][mp] + CHARGE_CURRENT_DN_DN[m][mp] + CHARGE_CURRENT_UP_DN[m][mp] + CHARGE_CURRENT_DN_UP[m][mp];
				Total_Spin_Current_at_site += CHARGE_CURRENT_UP_UP[m][mp] - CHARGE_CURRENT_DN_DN[m][mp] - CHARGE_CURRENT_UP_DN[m][mp] + CHARGE_CURRENT_DN_UP[m][mp];
				}
				}
			}
		}
		CC_Conservation_out<<m<<"    "<<Total_Charge_Current_at_site.real()<<"		"<<Total_Spin_Current_at_site.real()<<endl;
		}
	}
}


}
//-----------------------------------------------------------------------------------------//


//------------------------------Optical Conductivity---------------------------------------//
if(OPT_COND==true){
Mat_2_Complex_doub J1_MAT,J2_MAT;
int s_size=2,a_size=3;

J1_MAT.resize(H_size);
J2_MAT.resize(H_size);
for (int i=0;i<H_size;i++){
	J1_MAT[i].resize(H_size);
	J2_MAT[i].resize(H_size);
}


int x1,x2,temp_x1,temp_x2;
complex<double> iota(0.0,1.0);
complex<double> Zero_Complex(0.0,0.0);

for(int n=0;n<H_size;n++){
for(int m=0;m<H_size;m++){
	J1_MAT[n][m]=Zero_Complex;
	J2_MAT[n][m]=Zero_Complex;
}
}

complex<double> TEMP_J1_hop,TEMP_J2_hop,TEMP_J1_soc,TEMP_J2_soc;
int j,i1,i2,i3,i4,i5,i6,i2_p,i3_p,i5_p,i6_p;

for(int n=0;n<H_size;n++){
for(int m=0;m<H_size;m++){
TEMP_J1_hop=Zero_Complex;TEMP_J2_hop=Zero_Complex;
TEMP_J1_soc=Zero_Complex;TEMP_J2_soc=Zero_Complex;

for(int iy=0;iy<M;iy++){
for(int ix=0;ix<N;ix++){
for(int a=0;a<a_size;a++){

if(a==1){
j  = a   + 3*ix     + 3*N*iy;
i1 = a-1 + 3*ix     + 3*N*iy;
i4 = a+1 + 3*ix     + 3*N*iy;
i3 = a-1 + 3*(ix+1) + 3*N*iy;
i6 = a+1 + 3*(ix-1) + 3*N*iy;
i2 = a+1 + 3*ix + 3*N*(iy-1);
i5 = a-1 + 3*ix + 3*N*(iy+1);

i6_p = a+1 + 3*(N-1) + 3*N*iy;
i3_p = a-1 + 3*N*iy;
i2_p = a+1 + 3*ix + 3*N*(M-1);
i5_p = a-1 + 3*ix;

if(ix!=0 && ix!=N-1){
	TEMP_J1_hop += -iota*( (conj(Eigenvectors[n][j]))*(Eigenvectors[m][i1]) - (conj(Eigenvectors[n][i1]))*(Eigenvectors[m][j]) ) -iota*( (conj(Eigenvectors[n][j+MHS]))*(Eigenvectors[m][i1+MHS]) - (conj(Eigenvectors[n][i1+MHS]))*(Eigenvectors[m][j+MHS]) ) 
		   -iota*( (conj(Eigenvectors[n][j]))*(Eigenvectors[m][i6]) - (conj(Eigenvectors[n][i6]))*(Eigenvectors[m][j]) ) -iota*( (conj(Eigenvectors[n][j+MHS]))*(Eigenvectors[m][i6+MHS]) - (conj(Eigenvectors[n][i6+MHS]))*(Eigenvectors[m][j+MHS]) )
		   -iota*( (conj(Eigenvectors[n][i4]))*(Eigenvectors[m][j]) - (conj(Eigenvectors[n][j]))*(Eigenvectors[m][i4]) ) -iota*( (conj(Eigenvectors[n][i4+MHS]))*(Eigenvectors[m][j+MHS]) - (conj(Eigenvectors[n][j+MHS]))*(Eigenvectors[m][i4+MHS]) )
		   -iota*( (conj(Eigenvectors[n][i3]))*(Eigenvectors[m][j]) - (conj(Eigenvectors[n][j]))*(Eigenvectors[m][i3]) ) -iota*( (conj(Eigenvectors[n][i3+MHS]))*(Eigenvectors[m][j+MHS]) - (conj(Eigenvectors[n][j+MHS]))*(Eigenvectors[m][i3+MHS]) );
	
	TEMP_J1_soc += -LAMBDA_SOC*( (conj(Eigenvectors[n][j]))*(Eigenvectors[m][i1+MHS]) + (conj(Eigenvectors[n][i1+MHS]))*(Eigenvectors[m][j]) ) - LAMBDA_SOC*( (conj(Eigenvectors[n][j+MHS]))*(Eigenvectors[m][i1]) + (conj(Eigenvectors[n][i1]))*(Eigenvectors[m][j+MHS]) )
		       -LAMBDA_SOC*( exp(-2.0*PI*iota/3.0)*(conj(Eigenvectors[n][j]))*(Eigenvectors[m][i6+MHS]) + exp(2.0*PI*iota/3.0)*(conj(Eigenvectors[n][i6+MHS]))*(Eigenvectors[m][j]) ) - LAMBDA_SOC*( exp(2.0*PI*iota/3.0)*(conj(Eigenvectors[n][j+MHS]))*(Eigenvectors[m][i6]) + exp(-2.0*PI*iota/3.0)*(conj(Eigenvectors[n][i6]))*(Eigenvectors[m][j+MHS]) )
		       +LAMBDA_SOC*( (conj(Eigenvectors[n][i4]))*(Eigenvectors[m][j+MHS]) + (conj(Eigenvectors[n][j+MHS]))*(Eigenvectors[m][i4]) ) + LAMBDA_SOC*( (conj(Eigenvectors[n][i4+MHS]))*(Eigenvectors[m][j]) + (conj(Eigenvectors[n][j]))*(Eigenvectors[m][i4+MHS]) )
		       +LAMBDA_SOC*( exp(-2.0*PI*iota/3.0)*(conj(Eigenvectors[n][i3]))*(Eigenvectors[m][j+MHS]) + exp(2.0*PI*iota/3.0)*(conj(Eigenvectors[n][j+MHS]))*(Eigenvectors[m][i3]) ) + LAMBDA_SOC*( exp(2.0*PI*iota/3.0)*(conj(Eigenvectors[n][i3+MHS]))*(Eigenvectors[m][j]) + exp(-2.0*PI*iota/3.0)*(conj(Eigenvectors[n][j]))*(Eigenvectors[m][i3+MHS]) );
}

if(iy!=0 && iy!=M-1){
	TEMP_J2_hop += -iota*( (conj(Eigenvectors[n][j]))*(Eigenvectors[m][i1]) - (conj(Eigenvectors[n][i1]))*(Eigenvectors[m][j]) ) -iota*( (conj(Eigenvectors[n][j+MHS]))*(Eigenvectors[m][i1+MHS]) - (conj(Eigenvectors[n][i1+MHS]))*(Eigenvectors[m][j+MHS]) )
                   -iota*( (conj(Eigenvectors[n][j]))*(Eigenvectors[m][i2]) - (conj(Eigenvectors[n][i2]))*(Eigenvectors[m][j]) ) -iota*( (conj(Eigenvectors[n][j+MHS]))*(Eigenvectors[m][i2+MHS]) - (conj(Eigenvectors[n][i2+MHS]))*(Eigenvectors[m][j+MHS]) )
                   -iota*( (conj(Eigenvectors[n][i4]))*(Eigenvectors[m][j]) - (conj(Eigenvectors[n][j]))*(Eigenvectors[m][i4]) ) -iota*( (conj(Eigenvectors[n][i4+MHS]))*(Eigenvectors[m][j+MHS]) - (conj(Eigenvectors[n][j+MHS]))*(Eigenvectors[m][i4+MHS]) )
                   -iota*( (conj(Eigenvectors[n][i5]))*(Eigenvectors[m][j]) - (conj(Eigenvectors[n][j]))*(Eigenvectors[m][i5]) ) -iota*( (conj(Eigenvectors[n][i5+MHS]))*(Eigenvectors[m][j+MHS]) - (conj(Eigenvectors[n][j+MHS]))*(Eigenvectors[m][i5+MHS]) );

        TEMP_J2_soc += -LAMBDA_SOC*( (conj(Eigenvectors[n][j]))*(Eigenvectors[m][i1+MHS]) + (conj(Eigenvectors[n][i1+MHS]))*(Eigenvectors[m][j]) ) - LAMBDA_SOC*( (conj(Eigenvectors[n][j+MHS]))*(Eigenvectors[m][i1]) + (conj(Eigenvectors[n][i1]))*(Eigenvectors[m][j+MHS]) )
		       -LAMBDA_SOC*( exp(2.0*PI*iota/3.0)*(conj(Eigenvectors[n][j]))*(Eigenvectors[m][i2+MHS]) + exp(-2.0*PI*iota/3.0)*(conj(Eigenvectors[n][i2+MHS]))*(Eigenvectors[m][j]) ) - LAMBDA_SOC*( exp(-2.0*PI*iota/3.0)*(conj(Eigenvectors[n][j+MHS]))*(Eigenvectors[m][i2]) + exp(2.0*PI*iota/3.0)*(conj(Eigenvectors[n][i2]))*(Eigenvectors[m][j+MHS]) )
		       +LAMBDA_SOC*( (conj(Eigenvectors[n][i4]))*(Eigenvectors[m][j+MHS]) + (conj(Eigenvectors[n][j+MHS]))*(Eigenvectors[m][i4]) ) + LAMBDA_SOC*( (conj(Eigenvectors[n][i4+MHS]))*(Eigenvectors[m][j]) + (conj(Eigenvectors[n][j]))*(Eigenvectors[m][i4+MHS]) )
		       +LAMBDA_SOC*( exp(2.0*PI*iota/3.0)*(conj(Eigenvectors[n][i5]))*(Eigenvectors[m][j+MHS]) + exp(-2.0*PI*iota/3.0)*(conj(Eigenvectors[n][j+MHS]))*(Eigenvectors[m][i5]) ) + LAMBDA_SOC*( exp(-2.0*PI*iota/3.0)*(conj(Eigenvectors[n][i5+MHS]))*(Eigenvectors[m][j]) + exp(2.0*PI*iota/3.0)*(conj(Eigenvectors[n][j]))*(Eigenvectors[m][i5+MHS]) );
}

if(ix==0){
	TEMP_J1_hop += -iota*( (conj(Eigenvectors[n][j]))*(Eigenvectors[m][i1]) - (conj(Eigenvectors[n][i1]))*(Eigenvectors[m][j]) ) -iota*( (conj(Eigenvectors[n][j+MHS]))*(Eigenvectors[m][i1+MHS]) - (conj(Eigenvectors[n][i1+MHS]))*(Eigenvectors[m][j+MHS]) )
                   -iota*( (conj(Eigenvectors[n][i4]))*(Eigenvectors[m][j]) - (conj(Eigenvectors[n][j]))*(Eigenvectors[m][i4]) ) -iota*( (conj(Eigenvectors[n][i4+MHS]))*(Eigenvectors[m][j+MHS]) - (conj(Eigenvectors[n][j+MHS]))*(Eigenvectors[m][i4+MHS]) )
                   -iota*( (conj(Eigenvectors[n][i3]))*(Eigenvectors[m][j]) - (conj(Eigenvectors[n][j]))*(Eigenvectors[m][i3]) ) -iota*( (conj(Eigenvectors[n][i3+MHS]))*(Eigenvectors[m][j+MHS]) - (conj(Eigenvectors[n][j+MHS]))*(Eigenvectors[m][i3+MHS]) );

        TEMP_J1_soc += -LAMBDA_SOC*( (conj(Eigenvectors[n][j]))*(Eigenvectors[m][i1+MHS]) + (conj(Eigenvectors[n][i1+MHS]))*(Eigenvectors[m][j]) ) - LAMBDA_SOC*( (conj(Eigenvectors[n][j+MHS]))*(Eigenvectors[m][i1]) + (conj(Eigenvectors[n][i1]))*(Eigenvectors[m][j+MHS]) )
                       +LAMBDA_SOC*( (conj(Eigenvectors[n][i4]))*(Eigenvectors[m][j+MHS]) + (conj(Eigenvectors[n][j+MHS]))*(Eigenvectors[m][i4]) ) + LAMBDA_SOC*( (conj(Eigenvectors[n][i4+MHS]))*(Eigenvectors[m][j]) + (conj(Eigenvectors[n][j]))*(Eigenvectors[m][i4+MHS]) )
                       +LAMBDA_SOC*( exp(-2.0*PI*iota/3.0)*(conj(Eigenvectors[n][i3]))*(Eigenvectors[m][j+MHS]) + exp(2.0*PI*iota/3.0)*(conj(Eigenvectors[n][j+MHS]))*(Eigenvectors[m][i3]) ) + LAMBDA_SOC*( exp(2.0*PI*iota/3.0)*(conj(Eigenvectors[n][i3+MHS]))*(Eigenvectors[m][j]) + exp(-2.0*PI*iota/3.0)*(conj(Eigenvectors[n][j]))*(Eigenvectors[m][i3+MHS]) );
}

if(ix==N-1){
	TEMP_J1_hop += -iota*( (conj(Eigenvectors[n][j]))*(Eigenvectors[m][i1]) - (conj(Eigenvectors[n][i1]))*(Eigenvectors[m][j]) ) -iota*( (conj(Eigenvectors[n][j+MHS]))*(Eigenvectors[m][i1+MHS]) - (conj(Eigenvectors[n][i1+MHS]))*(Eigenvectors[m][j+MHS]) )
                   -iota*( (conj(Eigenvectors[n][j]))*(Eigenvectors[m][i6]) - (conj(Eigenvectors[n][i6]))*(Eigenvectors[m][j]) ) -iota*( (conj(Eigenvectors[n][j+MHS]))*(Eigenvectors[m][i6+MHS]) - (conj(Eigenvectors[n][i6+MHS]))*(Eigenvectors[m][j+MHS]) )
                   -iota*( (conj(Eigenvectors[n][i4]))*(Eigenvectors[m][j]) - (conj(Eigenvectors[n][j]))*(Eigenvectors[m][i4]) ) -iota*( (conj(Eigenvectors[n][i4+MHS]))*(Eigenvectors[m][j+MHS]) - (conj(Eigenvectors[n][j+MHS]))*(Eigenvectors[m][i4+MHS]) );

        TEMP_J1_soc += -LAMBDA_SOC*( (conj(Eigenvectors[n][j]))*(Eigenvectors[m][i1+MHS]) + (conj(Eigenvectors[n][i1+MHS]))*(Eigenvectors[m][j]) ) - LAMBDA_SOC*( (conj(Eigenvectors[n][j+MHS]))*(Eigenvectors[m][i1]) + (conj(Eigenvectors[n][i1]))*(Eigenvectors[m][j+MHS]) )
                       -LAMBDA_SOC*( exp(-2.0*PI*iota/3.0)*(conj(Eigenvectors[n][j]))*(Eigenvectors[m][i6+MHS]) + exp(2.0*PI*iota/3.0)*(conj(Eigenvectors[n][i6+MHS]))*(Eigenvectors[m][j]) ) - LAMBDA_SOC*( exp(2.0*PI*iota/3.0)*(conj(Eigenvectors[n][j+MHS]))*(Eigenvectors[m][i6]) + exp(-2.0*PI*iota/3.0)*(conj(Eigenvectors[n][i6]))*(Eigenvectors[m][j+MHS]) )
                       +LAMBDA_SOC*( (conj(Eigenvectors[n][i4]))*(Eigenvectors[m][j+MHS]) + (conj(Eigenvectors[n][j+MHS]))*(Eigenvectors[m][i4]) ) + LAMBDA_SOC*( (conj(Eigenvectors[n][i4+MHS]))*(Eigenvectors[m][j]) + (conj(Eigenvectors[n][j]))*(Eigenvectors[m][i4+MHS]) );
}

if(iy==0){
	TEMP_J2_hop += -iota*( (conj(Eigenvectors[n][j]))*(Eigenvectors[m][i1]) - (conj(Eigenvectors[n][i1]))*(Eigenvectors[m][j]) ) -iota*( (conj(Eigenvectors[n][j+MHS]))*(Eigenvectors[m][i1+MHS]) - (conj(Eigenvectors[n][i1+MHS]))*(Eigenvectors[m][j+MHS]) )
                   -iota*( (conj(Eigenvectors[n][i4]))*(Eigenvectors[m][j]) - (conj(Eigenvectors[n][j]))*(Eigenvectors[m][i4]) ) -iota*( (conj(Eigenvectors[n][i4+MHS]))*(Eigenvectors[m][j+MHS]) - (conj(Eigenvectors[n][j+MHS]))*(Eigenvectors[m][i4+MHS]) )
                   -iota*( (conj(Eigenvectors[n][i5]))*(Eigenvectors[m][j]) - (conj(Eigenvectors[n][j]))*(Eigenvectors[m][i5]) ) -iota*( (conj(Eigenvectors[n][i5+MHS]))*(Eigenvectors[m][j+MHS]) - (conj(Eigenvectors[n][j+MHS]))*(Eigenvectors[m][i5+MHS]) );

        TEMP_J2_soc += -LAMBDA_SOC*( (conj(Eigenvectors[n][j]))*(Eigenvectors[m][i1+MHS]) + (conj(Eigenvectors[n][i1+MHS]))*(Eigenvectors[m][j]) ) - LAMBDA_SOC*( (conj(Eigenvectors[n][j+MHS]))*(Eigenvectors[m][i1]) + (conj(Eigenvectors[n][i1]))*(Eigenvectors[m][j+MHS]) )
                       +LAMBDA_SOC*( (conj(Eigenvectors[n][i4]))*(Eigenvectors[m][j+MHS]) + (conj(Eigenvectors[n][j+MHS]))*(Eigenvectors[m][i4]) ) + LAMBDA_SOC*( (conj(Eigenvectors[n][i4+MHS]))*(Eigenvectors[m][j]) + (conj(Eigenvectors[n][j]))*(Eigenvectors[m][i4+MHS]) )
                       +LAMBDA_SOC*( exp(2.0*PI*iota/3.0)*(conj(Eigenvectors[n][i5]))*(Eigenvectors[m][j+MHS]) + exp(-2.0*PI*iota/3.0)*(conj(Eigenvectors[n][j+MHS]))*(Eigenvectors[m][i5]) ) + LAMBDA_SOC*( exp(-2.0*PI*iota/3.0)*(conj(Eigenvectors[n][i5+MHS]))*(Eigenvectors[m][j]) + exp(2.0*PI*iota/3.0)*(conj(Eigenvectors[n][j]))*(Eigenvectors[m][i5+MHS]) );
}

if(iy==M-1){
TEMP_J2_hop += -iota*( (conj(Eigenvectors[n][j]))*(Eigenvectors[m][i1]) - (conj(Eigenvectors[n][i1]))*(Eigenvectors[m][j]) ) -iota*( (conj(Eigenvectors[n][j+MHS]))*(Eigenvectors[m][i1+MHS]) - (conj(Eigenvectors[n][i1+MHS]))*(Eigenvectors[m][j+MHS]) )
                   -iota*( (conj(Eigenvectors[n][j]))*(Eigenvectors[m][i2]) - (conj(Eigenvectors[n][i2]))*(Eigenvectors[m][j]) ) -iota*( (conj(Eigenvectors[n][j+MHS]))*(Eigenvectors[m][i2+MHS]) - (conj(Eigenvectors[n][i2+MHS]))*(Eigenvectors[m][j+MHS]) )
                   -iota*( (conj(Eigenvectors[n][i4]))*(Eigenvectors[m][j]) - (conj(Eigenvectors[n][j]))*(Eigenvectors[m][i4]) ) -iota*( (conj(Eigenvectors[n][i4+MHS]))*(Eigenvectors[m][j+MHS]) - (conj(Eigenvectors[n][j+MHS]))*(Eigenvectors[m][i4+MHS]) );

        TEMP_J2_soc += -LAMBDA_SOC*( (conj(Eigenvectors[n][j]))*(Eigenvectors[m][i1+MHS]) + (conj(Eigenvectors[n][i1+MHS]))*(Eigenvectors[m][j]) ) - LAMBDA_SOC*( (conj(Eigenvectors[n][j+MHS]))*(Eigenvectors[m][i1]) + (conj(Eigenvectors[n][i1]))*(Eigenvectors[m][j+MHS]) )
                       -LAMBDA_SOC*( exp(2.0*PI*iota/3.0)*(conj(Eigenvectors[n][j]))*(Eigenvectors[m][i2+MHS]) + exp(-2.0*PI*iota/3.0)*(conj(Eigenvectors[n][i2+MHS]))*(Eigenvectors[m][j]) ) - LAMBDA_SOC*( exp(-2.0*PI*iota/3.0)*(conj(Eigenvectors[n][j+MHS]))*(Eigenvectors[m][i2]) + exp(2.0*PI*iota/3.0)*(conj(Eigenvectors[n][i2]))*(Eigenvectors[m][j+MHS]) )
                       +LAMBDA_SOC*( (conj(Eigenvectors[n][i4]))*(Eigenvectors[m][j+MHS]) + (conj(Eigenvectors[n][j+MHS]))*(Eigenvectors[m][i4]) ) + LAMBDA_SOC*( (conj(Eigenvectors[n][i4+MHS]))*(Eigenvectors[m][j]) + (conj(Eigenvectors[n][j]))*(Eigenvectors[m][i4+MHS]) );
}

}
/*
for(int j1=0;j1<M;j1++){
for(int i1=0;i1<N;i1++){
for(int a1=0;a1<a_size;a1++){
x1 = a1 + 3*i1 + 3*N*j1;

	for(int j2=0;j2<M;j2++){
	for(int i2=0;i2<N;i2++){
	for(int a2=0;a2<a_size;a2++){
	x2 = a2 + 3*i2 + 3*N*j2;
	
	if(i2==i1 && j2==j1){
		if(a1!=2 && a2==a1+1){
		TEMP_J1 += -iota*( (conj(Eigenvectors[n][x2]))*(Eigenvectors[m][x1]) - (conj(Eigenvectors[n][x1]))*(Eigenvectors[m][x2]) ) -iota*( (conj(Eigenvectors[n][x2+MHS]))*(Eigenvectors[m][x1+MHS]) - (conj(Eigenvectors[n][x1+MHS]))*(Eigenvectors[m][x2+MHS]) ) ;
		TEMP_J2 += -iota*( (conj(Eigenvectors[n][x2]))*(Eigenvectors[m][x1]) - (conj(Eigenvectors[n][x1]))*(Eigenvectors[m][x2]) ) -iota*( (conj(Eigenvectors[n][x2+MHS]))*(Eigenvectors[m][x1+MHS]) - (conj(Eigenvectors[n][x1+MHS]))*(Eigenvectors[m][x2+MHS]) ) ;
		}
		
		if(a1==0 && a2==a1+1){
		TEMP_J1 += -LAMBDA_SOC*( (conj(Eigenvectors[n][x2]))*(Eigenvectors[m][x1+MHS]) + (conj(Eigenvectors[n][x1+MHS]))*(Eigenvectors[m][x2]) ) - LAMBDA_SOC*( (conj(Eigenvectors[n][x2+MHS]))*(Eigenvectors[m][x1]) + (conj(Eigenvectors[n][x1]))*(Eigenvectors[m][x2+MHS]) );
		TEMP_J2 += -LAMBDA_SOC*( (conj(Eigenvectors[n][x2]))*(Eigenvectors[m][x1+MHS]) + (conj(Eigenvectors[n][x1+MHS]))*(Eigenvectors[m][x2]) ) - LAMBDA_SOC*( (conj(Eigenvectors[n][x2+MHS]))*(Eigenvectors[m][x1]) + (conj(Eigenvectors[n][x1]))*(Eigenvectors[m][x2+MHS]) );
		}

		if(a1==1 && a2==a1+1){
		TEMP_J1 += LAMBDA_SOC*( (conj(Eigenvectors[n][x2]))*(Eigenvectors[m][x1+MHS]) + (conj(Eigenvectors[n][x1+MHS]))*(Eigenvectors[m][x2]) ) + LAMBDA_SOC*( (conj(Eigenvectors[n][x2+MHS]))*(Eigenvectors[m][x1]) + (conj(Eigenvectors[n][x1]))*(Eigenvectors[m][x2+MHS]) );
		TEMP_J2 += LAMBDA_SOC*( (conj(Eigenvectors[n][x2]))*(Eigenvectors[m][x1+MHS]) + (conj(Eigenvectors[n][x1+MHS]))*(Eigenvectors[m][x2]) ) + LAMBDA_SOC*( (conj(Eigenvectors[n][x2+MHS]))*(Eigenvectors[m][x1]) + (conj(Eigenvectors[n][x1]))*(Eigenvectors[m][x2+MHS]) );
		}
	}

	if(i1!=N-1){
	if(i2==i1+1 && j2==j1){
		if(a1!=0 && a2==a1-1){
		TEMP_J1 += -iota*( (conj(Eigenvectors[n][x2]))*(Eigenvectors[m][x1]) - (conj(Eigenvectors[n][x1]))*(Eigenvectors[m][x2]) ) -iota*( (conj(Eigenvectors[n][x2+MHS]))*(Eigenvectors[m][x1+MHS]) - (conj(Eigenvectors[n][x1+MHS]))*(Eigenvectors[m][x2+MHS]) ) ;
		}

		if(a1==1 && a2==a1-1){
		TEMP_J1 += LAMBDA_SOC*( exp(-2.0*PI*iota/3.0)*(conj(Eigenvectors[n][x2]))*(Eigenvectors[m][x1+MHS]) + exp(2.0*PI*iota/3.0)*(conj(Eigenvectors[n][x1+MHS]))*(Eigenvectors[m][x2]) ) + LAMBDA_SOC*( exp(2.0*PI*iota/3.0)*(conj(Eigenvectors[n][x2+MHS]))*(Eigenvectors[m][x1]) + exp(-2.0*PI*iota/3.0)*(conj(Eigenvectors[n][x1]))*(Eigenvectors[m][x2+MHS]) );
		}

		if(a1==2 && a2==a1-1){
		TEMP_J1 += -LAMBDA_SOC*( exp(-2.0*PI*iota/3.0)*(conj(Eigenvectors[n][x2]))*(Eigenvectors[m][x1+MHS]) + exp(2.0*PI*iota/3.0)*(conj(Eigenvectors[n][x1+MHS]))*(Eigenvectors[m][x2]) ) - LAMBDA_SOC*( exp(2.0*PI*iota/3.0)*(conj(Eigenvectors[n][x2+MHS]))*(Eigenvectors[m][x1]) + exp(-2.0*PI*iota/3.0)*(conj(Eigenvectors[n][x1]))*(Eigenvectors[m][x2+MHS]) );
		}
	}
	}

	if(j1!=M-1){
	if(j2==j1+1 && i2==i1){
		if(a1!=0 && a2==a1-1){
		TEMP_J2 += -iota*( (conj(Eigenvectors[n][x2]))*(Eigenvectors[m][x1]) - (conj(Eigenvectors[n][x1]))*(Eigenvectors[m][x2]) ) -iota*( (conj(Eigenvectors[n][x2+MHS]))*(Eigenvectors[m][x1+MHS]) - (conj(Eigenvectors[n][x1+MHS]))*(Eigenvectors[m][x2+MHS]) ) ;
		}

		if(a1==1 && a2==a1-1){
		TEMP_J2 += LAMBDA_SOC*( exp(2.0*PI*iota/3.0)*(conj(Eigenvectors[n][x2]))*(Eigenvectors[m][x1+MHS]) + exp(-2.0*PI*iota/3.0)*(conj(Eigenvectors[n][x1+MHS]))*(Eigenvectors[m][x2]) ) + LAMBDA_SOC*( exp(-2.0*PI*iota/3.0)*(conj(Eigenvectors[n][x2+MHS]))*(Eigenvectors[m][x1]) + exp(2.0*PI*iota/3.0)*(conj(Eigenvectors[n][x1]))*(Eigenvectors[m][x2+MHS]) );
		}

		if(a1==2 && a2==a1-1){
		TEMP_J2 += -LAMBDA_SOC*( exp(2.0*PI*iota/3.0)*(conj(Eigenvectors[n][x2]))*(Eigenvectors[m][x1+MHS]) + exp(-2.0*PI*iota/3.0)*(conj(Eigenvectors[n][x1+MHS]))*(Eigenvectors[m][x2]) ) - LAMBDA_SOC*( exp(-2.0*PI*iota/3.0)*(conj(Eigenvectors[n][x2+MHS]))*(Eigenvectors[m][x1]) + exp(2.0*PI*iota/3.0)*(conj(Eigenvectors[n][x1]))*(Eigenvectors[m][x2+MHS]) );
		}
	}
	}

	}
	}
	}
*/	

}
}
}


J1_MAT[n][m]=TEMP_J1_hop + TEMP_J1_soc;
J2_MAT[n][m]=TEMP_J2_hop + TEMP_J2_soc;
}
}

double eta=0.02;
double mu_min=-5.0;
double mu_max=5.0;
double d_mu=0.005;
int mu_size = (int) ( (mu_max-mu_min)/d_mu );
double mu_;

complex<double> SIGMA_XY,SIGMA_XX,SIGMA_YY;


string OC_file="Optical_Conductivity_vs_mu.txt";
ofstream OC_file_out(OC_file.c_str());

for(int o=0;o<mu_size;o++){
        SIGMA_XY=Zero_Complex;
        SIGMA_XX=Zero_Complex;
        SIGMA_YY=Zero_Complex;
//      RHO_XY=Zero_Complex;
//        RHO_XX=Zero_Complex;

        mu_=mu_min + o*d_mu;

        for(int n=0;n<H_size;n++){
        for(int m=0;m<H_size;m++){

if(abs(EVALS[m]-EVALS[n])>1e-8){
        SIGMA_XY += ((2.0*PI)/((N*M*3.0)) )*( (Fermi_(EVALS[m],mu_)-Fermi_(EVALS[n],mu_)) )*( (1.0)/( (EVALS[m]-EVALS[n])*(EVALS[m]-EVALS[n]) + eta*eta)  )*J1_MAT[m][n]*J2_MAT[n][m];
        SIGMA_XX += ((2.0*PI*eta)/((N*M*3.0)) )*( (Fermi_(EVALS[m],mu_)-Fermi_(EVALS[n],mu_))/(EVALS[n]-EVALS[m]) )*( (1.0)/( (EVALS[m]-EVALS[n])*(EVALS[m]-EVALS[n]) + eta*eta)  )*J1_MAT[m][n]*J1_MAT[n][m];
        }
        }
        }

//RHO_XX = (SIGMA_YY)/( -(conj(SIGMA_XY)*SIGMA_XY) + SIGMA_XX*SIGMA_YY );
//RHO_XY = -(SIGMA_XY)/( -(conj(SIGMA_XY)*SIGMA_XY) + SIGMA_XX*SIGMA_YY );

OC_file_out<<mu_<<"     "<<SIGMA_XY.imag()<<"   "<<SIGMA_XX.real()<<endl;
}

}

/*
double eta=0.1;
double w_min=-5;
double w_max=5;
double dw=0.0025;
int w_size = (int) ( (w_max-w_min)/dw );
double omega;

complex<double> SIGMA_XY,SIGMA_XX,SIGMA_YY;
complex<double> RHO_XY,RHO_XX;
omega=w_min;

//string OC_file="Optical_Conductivity_vs_omega.txt";
//ofstream OC_file_out(OC_file.c_str());

//for(int o=0;o<w_size;o++){
	SIGMA_XY=Zero_Complex;
	SIGMA_XX=Zero_Complex;
	SIGMA_YY=Zero_Complex;
	RHO_XY=Zero_Complex;
        RHO_XX=Zero_Complex;
//	omega=omega + dw;

	for(int n=0;n<H_size;n++){
	for(int m=0;m<H_size;m++){

	if(EVALS[m]!=EVALS[n]){
//	SIGMA += ((1.0)/(sqrt(N*M*1.0)) )*( (Fermi_(EVALS[n],mu)-Fermi_(EVALS[m],mu))/(EVALS[m]-EVALS[n]) )*( (1.0*eta)/( (omega + EVALS[n]-EVALS[m])*(omega + EVALS[n]-EVALS[m]) + eta*eta )  )*J1_MAT[n][m]*J2_MAT[m][n];
	SIGMA_XY += ((2.0*PI)/((N*N*3.0)) )*( (Fermi_(EVALS[m],mu)-Fermi_(EVALS[n],mu)) )*( (1.0)/( (EVALS[m]-EVALS[n])*(EVALS[m]-EVALS[n]) + eta*eta)  )*J1_MAT[m][n]*J2_MAT[n][m];
	SIGMA_XX += ((2.0*PI*eta)/((N*N*1.0)) )*( (Fermi_(EVALS[m],mu)-Fermi_(EVALS[n],mu))/(EVALS[n]-EVALS[m]) )*( (1.0)/( (EVALS[m]-EVALS[n])*(EVALS[m]-EVALS[n]) + eta*eta)  )*J1_MAT[m][n]*J1_MAT[n][m];
	SIGMA_YY += ((2.0*PI*eta)/((M*M*1.0)) )*( (Fermi_(EVALS[m],mu)-Fermi_(EVALS[n],mu))/(EVALS[n]-EVALS[m]) )*( (1.0)/( (EVALS[m]-EVALS[n])*(EVALS[m]-EVALS[n]) + eta*eta)  )*J2_MAT[m][n]*J2_MAT[n][m];

	}
	}
	}

RHO_XX = (SIGMA_YY)/( -(conj(SIGMA_XY)*SIGMA_XY) + SIGMA_XX*SIGMA_YY );
RHO_XY = -(SIGMA_XY)/( -(conj(SIGMA_XY)*SIGMA_XY) + SIGMA_XX*SIGMA_YY );

cout<<EVALS[Np_-1]<<"	"<<mu<<"	"<<B_mag<<"	"<<SIGMA_XY.imag()<<"	"<<SIGMA_XX.real()<<"	"<<RHO_XX.real()<<"	"<<RHO_XY.imag()<<endl;

//OC_file_out<<omega<<"		"<<SIGMA.real()<<"	"<<SIGMA.imag()<<endl;

//}

}
*/

return 0;
}
