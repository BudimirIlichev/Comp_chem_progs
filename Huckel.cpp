#include <iostream>
#include <cmath>
#include <vector>
#include <numbers>
#include "Chem_Table.h"

using namespace std;

double  BondOrder(int i,int j, int n ,const vector<double>& Popul,const vector<vector<double>>& V);

int main(){
    int n,i,j,num_b;
    double phi,c,s;
    double pi=acos(-1.0);
    string symbol;
    cout<<"Number of atoms:\n";
    cin>>n; 
    cout<<"Number of  bond:\n";
    cin>>num_b;

    ChemicalTable table;

    vector<vector<double>> F(n,vector<double>(n));
    vector<vector<double>> V(n,vector<double>(n));
    vector<double> E(n,0.0);

    for(i=0;i<n;i++){
        while(true){
        cout<<"Name of element"<<endl;
        cin>>symbol;
        if (table.hasElement(symbol)){
            F[i][i]=table.getHuckelh(symbol);
            break;
        }else{
            cout<<"Element is not found"<<endl;
        }
    }
}
    for(int k=0;k<num_b;k++){
        while(true){
            int atom1,atom2;
            cout<<"Write number of connected atoms:"<<endl;
            cin>>atom1>>atom2;
            i=atom1-1;
            j=atom2-1;
            if(i>=0 && i<n && j>=0 && j<n && i!=j){
                F[i][j]=1.0;
                F[j][i]=1.0;
                break;
            }else{
                cout<<"Wrong indexes, repeat input"<<endl;
            }
        }
    }
    for(i=0;i<n;i++){
        V[i][i]=1.0;
    }
    double Eps=1e-9;
    while(true){

     double max_val=0.0;
     double p=0;
     double q=0;
    
    
    for (i=0;i<n-1;i++){
        for(j=i+1;j<n;j++){
            if (std::abs(F[i][j])>max_val){
                max_val=std::abs(F[i][j]);
                p=i;
                q=j;
            }
        }
    }

            if (max_val<Eps){
                break;
            }
            if ((std::abs(F[p][p]-F[q][q])<1e-15)){
                 phi=pi/4;
            }else{
                phi=0.5*atan2(2*F[p][q],(F[q][q]-F[p][p]));
            }
             c=cos(phi);
             s=sin(phi);

    for(i=0;i<n;i++){
        if(i!=p && i!=q){
            double F_ip_old=F[i][p];
            double F_iq_old=F[i][q];

            F[i][p]=F_ip_old*c-F_iq_old*s;
            F[i][q]=F_iq_old*c+F_ip_old*s;

            F[p][i]=F[i][p];
            F[q][i]=F[i][q];
        }
    }
    double F_pp_old=F[p][p];
    double F_qq_old=F[q][q];
    double F_pq_old=F[p][q];

    F[p][p]=F_pp_old*c*c-2*F_pq_old*c*s+F_qq_old*s*s;
    F[q][q]=F_qq_old*c*c+2*F_pq_old*c*s+F_pp_old*s*s;

    F[p][q]=0.0;
    F[q][p]=0.0;

    for(i=0;i<n;i++){
        double V_old_ip=V[i][p];
        V[i][p]=c*V_old_ip-s*V[i][q];
        V[i][q]=s*V_old_ip+c*V[i][q];
    }
    }
    for(i=0;i<n;i++){
        E[i]=F[i][i];
    }
    for(i=0;i<n;i++){
        for(j=0;j<n-i-1;j++){
            if(E[j]<E[j+1]){
                 std::swap(E[j],E[j+1]);
                  for(int k=0;k<n;k++){
                    std::swap(V[k][j],V[k][j+1]);
                  }
                }
            }   
    }
   vector<double>Popul(n,0.0);
   int rem_el=n;
   for(i=0;i<n;i++){
    if(rem_el>=2){
        Popul[i]=2;
        rem_el-=2;
    }else if (rem_el==1){
        Popul[i]=1;
        rem_el-=1;
    }
   }
   double Energy_beta=0.0;
   for(i=0;i<n;i++){
    Energy_beta+=Popul[i]*E[i];
   }
   double p_01=BondOrder(0,1,n,Popul,V);
   cout<<"Bond Order:" <<p_01<<endl;
}
double  BondOrder(int i,int j, int n ,const vector<double>& Popul,const vector<vector<double>>& V){
    double p_ij=0.0;
    for(int k=0;k<n;k++){
        p_ij+=Popul[k]*V[i][k]*V[j][k];
    }
    return p_ij;
}
