//
// Created by User on 2022-03-10.
//
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "FEM.h"
/*
Materials NewMaterials(){


};
*/

void read_input(){

    char charr[256];
    int i, j, temp1;
    ifstream fin;
    cout << "--------------------  input  --------------------" << endl;
    cout << endl;
    fin.open("input - 1d.txt");
    fin.getline(charr, 256);
    fin >> temp1;
    cout <<temp1;
    fin.close();
    cout << "--------------------  -----  --------------------" << endl;
    cout << endl;

}