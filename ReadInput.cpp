//
// Created by User on 2022-03-10.
//
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "FEM.h"
using namespace std;
/*
Materials NewMaterials(){


};
*/
// 고유빈바보
void read_input(){

    int count_line = 0;

    ifstream fin;
    string buffer;
    // print input
    cout << "--------------------  input  --------------------" << endl;
    cout << endl;

    fin.open("input-1d.txt");

    while (fin.peek() != EOF) {
        getline(fin, buffer);
        cout << buffer <<endl;
        count_line++; // total lines
    }
    fin.close();
    // data save
    fin.open("input-1d.txt");
    double a[count_line];
    for(int i=0; i<count_line; i++){
        fin>>a[i];
    }
    fin.close();
    for(int i=0; i < a[1];i++){

        Materials Mat[i];
        // Mat[i].Sig_t= a[];
    }
    for(int i =0;i <count_line;i++){
        cout<<to_string(a[i])<<endl;
    }
}
