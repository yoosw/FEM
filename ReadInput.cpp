//
// Created by User on 2022-03-10.
//
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "FEM.h"

Materials NewMaterials(){
    //int num, group;
    //    double Sig_t, Sig_a, Sig_ii, Sig_ij, Nusig_f;

}

double read_input(){
    char charr[256];
    int BCon[2];
    int i, j, k1, k2;

    ifstream fin;

    cout << "--------------------  input  --------------------" << endl;
    cout << endl;
    fin.open("input - 1d.txt");

    fin.getline(charr, 256);
    fin >> Materials;

    fin.getline(charr, 256);
    fin.getline(charr, 256);
    cout << charr << endl;
    fin >> Reg_X;
    Reg_Dis_X = new double[Reg_X];
    cout<<"Region division:" << Reg_X << endl;

    fin.getline(charr, 256);
    fin.getline(charr, 256);
    cout << charr << endl;
    for (i = 0; i < Reg_X; i++)
    {
        fin >> Reg_Dis_X[i];
    }

    for (i = 0; i < Reg_X; i++)
    {
        cout << Reg_Dis_X[i];
    }
    cout << endl;


    fin.getline(charr, 256);
    fin.getline(charr, 256);
    cout << charr << endl;
    fin >> node_length;
    cout << node_length << endl;
    Tot_X = Reg_Dis_X[Reg_X-1] / node_length;
    hx = 1 / Tot_X;
    Tot_Y = 1;
    Mat = new int[Reg_X];


    fin.getline(charr, 256);
    fin.getline(charr, 256);
    cout << charr << endl;

    for (i = 0; i < Reg_X; i++)
    {
        fin >> Mat[i];
        cout << Mat[i]<<" ";
    }
    cout<< endl;

    Sig_t = new double* [Materials];
    Sig_a = new double* [Materials];
    Sig_ii = new double* [Materials];
    Sig_ij = new double* [Materials];
    NuSig_f = new double* [Materials];

    for (i = 0; i < Materials; i++) {
        Sig_t[i] = new double[2];
        Sig_a[i] = new double[2];
        Sig_ii[i] = new double[2];
        Sig_ij[i] = new double[2];
        NuSig_f[i] = new double[2];
    }

    fin.getline(charr, 256);

    int count1 = 1;
    for (i = 0; i < Materials; i++) {
        for (j = 0; j < 2; j++)
        {
            if (count1 ==1) {
                fin.getline(charr, 256);
                cout << charr << endl;
                fin >> k1 >> k2 >> Sig_t[i][j] >> Sig_a[i][j] >> Sig_ii[i][j] >> Sig_ij[i][j] >> NuSig_f[i][j];
                cout << k1 << "\t" << k2 << "\t" << Sig_t[i][j] << "\t" << Sig_a[i][j] << "\t" << Sig_ii[i][j] << "\t" << Sig_ij[i][j] << "\t" << NuSig_f[i][j];
                cout << endl;
            }
            if (count1 > 1) {
                fin.getline(charr, 256);
                fin >> k1 >> k2 >> Sig_t[i][j] >> Sig_a[i][j] >> Sig_ii[i][j] >> Sig_ij[i][j] >> NuSig_f[i][j];
                cout << k1 << "\t" <<k2 << "\t" << Sig_t[i][j] << "\t" << Sig_a[i][j] << "\t" << Sig_ii[i][j] << "\t" << Sig_ij[i][j] << "\t" << NuSig_f[i][j];
                cout << endl;
            }
            count1 = count1 + 1;
        }
    }


    D1 = 1 / Sig_t[0][0];
    D2 = 1 / Sig_t[1][0];

    fin.getline(charr, 256);
    fin.getline(charr, 256);
    cout << charr << endl;

    fin >> BCon[0] >> BCon[1];
    cout << BCon[0] << ' ' << BCon[1] << ' ' << endl;

    fin.getline(charr, 256);
    fin.getline(charr, 256);
    cout << charr << endl;

    fin >> eps_k >> eps_phi;
    cout << eps_k << ' ' << eps_phi << endl;

    fin.close();
    cout << "--------------------  -----  --------------------" << endl;
    cout << endl;

}