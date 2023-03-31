/*
    adaptive_ising.cpp
    Supporting material software to the article 
    "Statistical modeling of adaptive neural networks 
    explains coexistence of avalanches and oscillations in resting human brain "
    Copyright (C) 2022 F. Lombardi et al.
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
                                                                           */

#include "iostream"
#include "stdlib.h"
#include "math.h"
#include "time.h"
#include "fstream"
#include "vector"
#include "string"

using namespace std;

#define N 10000   // size of the systems
int s[N];          // spin variables
double m;         // magnetization
double H;           //external field
double beta=0.987;              // inverse temperature
// double beta=0.98;              // inverse temperature
double dt = 1./double(N);          // time update 
double c = 0.01;              // feedback strength

double casual(){
// a random number uniform in (0,1)
    double   x = double  (rand())/(RAND_MAX+1.0);
    return x;
}
void init(){
// initialization: infinite temperature configuration and H=0
    for(int i=0;i<=N-1;i++){                    
        s[i]=1;
        if(i%2==0) s[i]=-1;
        m+=double(s[i])/double(N);
    }
    H=0;
}
void single_update(int n){
// update of the spin n according to the heat bath method
    int news;
    double pp = 1./(1.+1./exp(2.*beta*(m+H)));
    if(casual()<pp) news=1;
    else news =-1;
    if(news!=s[n]){
        s[n]=-s[n];
        m += double(2*s[n])/double(N);
    }
    H-=c*m*dt; // I update also the external field according to the feedback
}
void update(){
// one sweep over N spins
    for(int i=0;i<=N-1;i++) {
        int n =int(N*casual());
        single_update(n);
    }
}

int main (int argc, char **argv){
    beta = atof(argv[1]);
    srand(time(0));              // seed for the random number generator
    int t=0;                     // time 1=sweep over N spins
    int wait =100;               // waiting time before printing/measuring
    int tot =10000;           // total time of the simulation             1=one sweep over N spins
    init();
    do{
        if(t>wait)  cout << m << "  " << H << endl;          // print m and H
        update();
        t++;
       } while(t<=wait + tot);
    return 0 ;
}