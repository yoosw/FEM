//
// Created by User on 2022-03-10.
//

#ifndef FEM_FEM_H
#define FEM_FEM_H
/*
int Reg_X, Reg_Y, Reg_Z, Tot_X, Tot_Y, Tot_Z;
double* Reg_Dis_X, * K_old, * K;
double node_length;
double hx;
double D1,D2;
int* Reg_Sub_X;
struct error;
struct Matrix;
*/
struct Materials;
struct Materials{

    int num, group;
    double Sig_t, Sig_a, Sig_ii, Sig_ij, Nusig_f;
    int BCon;
} typedef Materials;
//Materials NewMaterials();

void read_input();
void loadCrossSection();


using namespace std;
#endif //FEM_FEM_H
