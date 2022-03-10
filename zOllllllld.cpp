#include <iostream>
#include <fstream>
#include <cmath>




int Materials, Reg_X, Reg_Y, Reg_Z, Tot_X, Tot_Y, Tot_Z;
double* Reg_Dis_X, * K_old, * K;
double node_length;
double hx;
double D1,D2;
int* Reg_Sub_X;
int* Mat;
double** Sig_t, ** Sig_a, ** Sig_ii, ** Sig_ij, ** NuSig_f;
double eps_k, eps_phi;
int Dimension;
using namespace std;
void read_input()
{
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
double func_a1(int i_num, int k) {
    return hx*(-D1*(4*pow(i_num,2)-8/3)+(Sig_a[k][0]+Sig_ij[k][0])*(2*pow(i_num,4)+16/15));
}
double func_a2(int i_num, int k) {
    return hx * (-D2 * (4 * pow(i_num, 2) - 8 / 3) + (Sig_a[k][1]) * (2 * pow(i_num, 4) + 16 / 15));
}
double func_b1(int i_num, int k) {
    return hx * (-D1 * (4 * pow(i_num, 2) - 8 / 3) + (Sig_a[k][0] + Sig_ij[k][0]) * (2 * pow(i_num, 4) + 4 * pow(i_num,3)+2*pow(i_num,2) - 4 / 15));
}
double func_b2(int i_num, int k) {
    return hx * (-D2 * (4 * pow(i_num, 2) - 8 / 3) + (Sig_a[k][1]) * (2 * pow(i_num, 4) + 4 * pow(i_num, 3) + 2 * pow(i_num, 2) - 4 / 15));
}
double func_c1(int i_num, int k) {
    return hx * (-D1 * (4 * pow(i_num, 2) + 8 * i_num + 4 / 3) + (Sig_a[k][0] + Sig_ij[k][0]) * (2 * pow(i_num, 4) + 4 * pow(i_num, 3) + 2 * pow(i_num, 2) - 4 / 15));
}
double func_c2(int i_num, int k) {
    return hx * (-D2 * (4 * pow(i_num, 2) + 8 * i_num+4 / 3) + (Sig_a[k][1]) * (2 * pow(i_num, 4) + 4 * pow(i_num, 3) + 2 * pow(i_num, 2) - 4 / 15));
}
double func_d1(int i_num, int k) {
    return hx * (-D1 * (4 * pow(i_num, 2) + 8 * i_num+4 / 3) + (Sig_a[k][0] + Sig_ij[k][0]) * (2 * pow(i_num, 4)+8*pow(i_num,3)+12*pow(i_num,2)+8*i_num + 46 / 15));
}
double func_d2(int i_num, int k) {
    return hx * (-D2 * (4 * pow(i_num, 2) + 8 * i_num+4 / 3) + (Sig_a[k][1]) * (2 * pow(i_num, 4) + 8 * pow(i_num, 3) + 12 * pow(i_num, 2) + 8 * i_num + 46 / 15));
}
double func_e1(int i_num) {
    return hx * (2 * pow(i_num, 2) - 4 / 3);
}
double func_e2(int i_num) {
    return hx * (2 * pow(i_num, 2) - 4 / 3);
}
double func_f1(int i_num) {
    return hx * (2 * pow(i_num, 2)+ 4 * i_num + 2/3) ;
}
double func_f2(int i_num) {
    return hx * (2 * pow(i_num, 2) + 4 * i_num + 2 / 3);
}
double func_a1_d1(int i_num, int k) {
    return 2* hx * (D1 * (4 * pow(i_num, 2) + 8 * i_num+ 4 / 3) + (Sig_a[k][0] + Sig_ij[k][0]) * (2 * pow(i_num, 4) + 8 * pow(i_num, 3) + 12 * pow(i_num, 2) + 8 * i_num+46 / 15));
}
double func_a2_d2(int i_num, int k) {
    return 2* hx * (D2 * (4 * pow(i_num, 2) + 8 * i_num+4 / 3) + (Sig_a[k][1]) * (2 * pow(i_num, 4) + 8 * pow(i_num, 3) + 12 * pow(i_num, 2) + 8 * i_num+46 / 15));
}
double func_e1_f1(int i_num) {
    return hx * (4 * pow(i_num, 2) + 4 * i_num - 2 / 3);
}
double func_e2_f2(int i_num) {
    return hx * (4 * pow(i_num, 2) + 4 * i_num - 2 / 3);
}
double** AA, * X, * B;
double** NA1, ** NA2, * NB1, * NB2, * NB1_old, * NB2_old, * NC1, * NC2;
/*
void Jacobian_Iteration(double** AA, double* XX, double* xo, double* BB, int n)
{
    int i, j, Iteration;E
    double a, b, er1;

    Iteration = 1;

    do {
        for (i = 0; i < n; i++)
            xo[i] = XX[i];

        for (i = 0; i < n; i++) {
            er1 = 0.;
            a = 0.;
            b = 0.;
            for (j = 0; j < i; j++) {
                //a += AA[i][j] * xo[j];
                a += AA[i][j] * XX[j];

            }
            for (j = i + 1; j < n; j++) {
                b += AA[i][j] * xo[j];
            }
            XX[i] = (BB[i] - a - b) / AA[i][i];

            if ((abs(XX[i] - xo[i] / XX[i]) > er1) && (XX[i] != 0.))
                er1 = abs((XX[i] - xo[i] / XX[i]));
        }
        cout << "iteration: " << Iteration << "Error= " << er1 << endl;
        Iteration += 1;
    } while (er1 < 0.0001);
}

*/
void Gauss_Seidel_Iteration(double** AA, double* XX, double* xo, double* BB, int n)
{
    int i, j, Iteration;
    double a, b, er1;

    Iteration = 1;
    for (i = 0; i < n; i++) {
        xo[i] = XX[i];
    }


    while (Iteration<=17){
        for (i = 0; i < n; i++) {
            er1 = 0.;
            a = 0.;
            b = 0.;
            for (j = 0; j < i-1; j++) {
                //a += AA[i][j] * xo[j];
                a += AA[i][j] * XX[j];

            }
            for (j = i; j < n; j++) {
                b += AA[i][j] * xo[j];
            }
            XX[i] = (BB[i] - a - b) / AA[i][i];
            er1 = fabs((XX[i] - xo[i] / XX[i]));
            if (er1 < 0.00001) {
                cout << "iteration: " << Iteration << "Error= " << er1 << endl;
                break;
            }
        }
        cout << "iteration: " << Iteration << "Error= " << er1 << endl;
        Iteration += 1;
        for (i = 0; i < n; i++) {
            xo[i] = XX[i];
        }

    }

}

void make_matrix()
{

    int i, j;
    NA1 = new double* [Tot_X];
    NA2 = new double* [Tot_X];
    NB1 = new double[Tot_X];
    NB2 = new double[Tot_X];
    NC1 = new double[Tot_X];
    NC2 = new double[Tot_X];
    NB1_old = new double[Tot_X];
    NB2_old = new double[Tot_X];
    K = new double[Tot_X];
    K_old = new double[Tot_X];

    for (i = 0; i < Tot_X; i++) {
        NA1[i] = new double[Tot_X];
        NA2[i] = new double[Tot_X];
    }

    for (i = 0; i < Tot_X; i++) {

        NB1[i] = 0.;
        NB2[i] = 0.;
        NC1[i] = 0.;
        NC2[i] = 0.;
        NB1_old[i] = 0.;
        NB2_old[i] = 0.;
        for (j = 0; j < Tot_X; j++)
        {
            NA1[i][j] = 0.;
            NA2[i][j] = 0.;
        }

    }


    for (i = 0; i < Tot_X; i++)
    {
        for (j = 0; j < Tot_Y; j++)
        {

            for (int t = 0; t < Reg_X; t++)
            {
                if ((i < Reg_Dis_X[t] / node_length) && (j == i))
                {
                    NA1[i][j] = func_a1_d1(i, Mat[t]);
                    NA2[i][j] = func_a2_d2(i, Mat[t]);

                }

                if ((i < Reg_Dis_X[t] / node_length) && (j == i + 1))
                {
                    NA1[i][j] = func_b1(i, Mat[t]);
                    NA2[i][j] = func_b2(i, Mat[t]);

                }
                if ((i > 0) && (i < Reg_Dis_X[t] / node_length) && (j == i - 1))
                {
                    NA1[i][j] = func_c1(i, Mat[t]);
                    NA2[i][j] = func_c2(i, Mat[t]);

                }

            }

            if ((i == 0) && (j == 0))
            {
                NA1[i][j] = func_a1(i, Mat[0]);
                NA2[i][j] = func_a2(i, Mat[0]);

            }
            if ((i == 0) && (j == 1))
            {
                NA1[i][j] = func_b1(i, Mat[0]);
                NA2[i][j] = func_b2(i, Mat[0]);

            }
            if ((i == Tot_X) && (j == Tot_X))
            {
                NA1[i][j] = func_d1(i, Mat[Reg_X - 1]);
                NA2[i][j] = func_d2(i, Mat[Reg_X - 1]);

            }
            if ((i == Tot_X) && (j == Tot_X - 1))
            {
                NA1[i][j] = func_c1(i, Mat[Reg_X - 1]);
                NA2[i][j] = func_c2(i, Mat[Reg_X - 1]);

            }
        }
    }
    for (i = 0; i < Tot_X; i++)
    {
        NC1[i] = func_e1_f1(i);
        NC2[i] = func_e2_f2(i);
        if (i == 0)
        {
            NC1[i] = func_e1(i);
            NC2[i] = func_e2(i);

        }
        if (i == Tot_X)
        {
            NC1[i] = func_f1(i);
            NC2[i] = func_f2(i);
        }
    }

}
double c1, c2;
void update_keff()
{
    c1 = 0;
    c2 = 0;

    int i, j, t;
    for (i = 0; i < Tot_X; i++) {
        K_old[i] = K[i];
        for (j = 0; j < Tot_Y; j++) {
            t = Mat[i] - 1;
            if (t == 0) {
                c1 = c1 + (NuSig_f[0][t] * NB1[i] + NuSig_f[1][t] * NB2[i]);
                c2 = c2 + (NuSig_f[0][t] * NB1_old[i] + NuSig_f[1][t] * NB2_old[i]);
            }
            if (t == 1) {
                c1 = c1;
                c2 = c2;
            }
        }

        K[i] = float(K_old[i] * (c1 / c2));
    }

}
void print_out()
{
    for (int i = 0; i < Tot_X; i++) {
        cout << K[i] << endl;
    }

}
int main() {

    double kerr = 1., ferr = 1.;
    for (int i = 0; i < Tot_X; i++) {
        K[i] = 1;
        K_old[i] = 1;
    }

    read_input();
    make_matrix();
    /*
    do {
        Gauss_Seidel_Iteration(NA1, NB1, NB1_old, NC1, Tot_X);
        for (int i = 0; i < Tot_X; i++) {
            for (int j = 0; j < Tot_X; j++) {
                cout << NA1[i][j] << " ";
            }
            cout << endl;

        }
        Gauss_Seidel_Iteration(NA2, NB2, NB2_old, NC2, Tot_X);
        update_keff();

    } while ((kerr > eps_k) || (ferr > eps_phi));
    */
    for (int i = 0; i < Tot_X; i++) {
        for (int j = 0; j < Tot_X; j++) {
            cout << NA1[i][j] << " ";
        }
        cout << endl;

    }
    Gauss_Seidel_Iteration(NA1, NB1, NB1_old, NC1, Tot_X);
    Gauss_Seidel_Iteration(NA2, NB2, NB2_old, NC2, Tot_X);
    update_keff();
    print_out();
    return 0;
}
