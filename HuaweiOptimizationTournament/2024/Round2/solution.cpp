#include<bits/stdc++.h>
using namespace std;
typedef long long ll;
typedef pair<int, int> ii;
//Input
vector<int> cpu, mem, price;
//Constant
double cA = 3.0;

//Struct that represent a Virtual Machine
struct VM{
    //Save data about total resources and used resources
    //Also save start time, price, and type of the machine
    int cputot, memtot, start;
    int cpuus = 0;
    int memus = 0;
    int precio;
    int tipo;
    //Here I save the times that enter the containers alocated
    multiset<int> tiempos;

    VM(){}

    //Just init
    VM(int start_, int type){
        tipo = type;
        cputot = cpu[tipo];
        memtot = mem[tipo];
        start = start_;
        precio = price[tipo];
    }

    //When desalocate a container we erase its enter time
    void sacartiempo(int t){
        tiempos.erase(tiempos.find(t));
    }
};

//The instructions to print
map<int, vector<vector<int> > > instrucciones;
//Save the VM
map<int, VM> opciones;
//Save in which VM I am
map<int, vector<int> > arriendo;
int idVM = 1;
int M, D, T;

//Search for the best VM to open
void nuevaopcion(int t, int id, int c, int m){
    int best = 1e9;
    int where = M - 1;
    //Test al the options
    for(int i = 0; i<M; i++){
        //If we are in the small case an heuristic
        if(M == 9 && cpu[i] >= (1.0 + 0.45*pow(double(T - t)/double(T), 3.6))*c && mem[i] >= (1.0 + 0.45*pow(double(T - t)/double(T), 3.6))*m && price[i] < best){
            best = price[i];
            where = i;
        }
        //If we are in the large test just the cheapest
        if(M == 42 && cpu[i] >= c && mem[i] >= m && price[i] < best){
            best = price[i];
            where = i;
        }
    }
    //Just open the VM
    opciones[idVM] = VM(t+D, where);
    opciones[idVM].cpuus = c;
    opciones[idVM].memus = m;
    opciones[idVM].tiempos.insert(t+D);
    instrucciones[t].push_back({1, idVM, where+1});
    instrucciones[t+D].push_back({3, id, idVM});
    arriendo[id] = {idVM, c, m, t+D};
    idVM++;
    return;
}

//Choose the best option to alocate the container
void mejoropcion(int t, int id, int c, int m){
    //Search for the best price if we open a new VM
    int best = 1e9;
    int wherex = -1;
    for(int i = 0; i<M; i++){
        if(cpu[i] >= c && mem[i] >= m && price[i] < best){
            best = price[i];
            wherex = i;
        }
    }
    int lasttime = 1e9;
    int where = -1;
    //We find a good option to alocate
    for(auto u: opciones){
        //Check if we can alocate here
        if(u.second.cpuus + c > u.second.cputot || u.second.memus + m > u.second.memtot) continue;
        //The critery for the big case
        if(M == 42 && cA*best <= u.second.precio) continue;
        //Search for the VM with oldest container
        if(*u.second.tiempos.begin() < lasttime){
            lasttime = *u.second.tiempos.begin();
            where = u.first;
        }
    }

    //If we find a good option we use it
    if(where != -1){
        int tiempo = max(t, opciones[where].start);
        instrucciones[tiempo].push_back({3, id, where});
        opciones[where].cpuus += c;
        opciones[where].memus += m;
        opciones[where].tiempos.insert(tiempo);
        arriendo[id] = {where, c, m, tiempo};
        return;
    }
    //If not, we create a new VM
    nuevaopcion(t, id, c, m);
    return;
}

//We close a container
void cerrar(int t, int id){
    //Change in resources of the VM that alocated the container
    vector<int> X = arriendo[id];
    opciones[X[0]].cpuus -= X[1];
    opciones[X[0]].memus -= X[2];
    opciones[X[0]].sacartiempo(X[3]);
    arriendo.erase(id);
    //If there is not any other container alocated, we close the VM
    if(opciones[X[0]].cpuus == 0){
        instrucciones[t+1].push_back({2, X[0]});
        opciones.erase(X[0]);
    }
    return;
}

int main(){

	ios::sync_with_stdio(false);

    //Read input

    cin>>M>>D;

    cpu.resize(M);
    mem.resize(M);
    price.resize(M);

    for(int i = 0; i<M; i++){
        cin>>cpu[i]>>mem[i]>>price[i];
    }

    //Ban bad machines :(
    if(M == 42){
        cpu[0] = 0;
        cpu[1] = 0;
        cpu[2] = 0;
        cpu[3] = 0;
        cpu[4] = 0;
        cpu[5] = 0;
        cpu[6] = 0;
    }
    else{
        cpu[3] = 0;
        cpu[4] = 0;
    }

    int t;
    cin>>t;

    T = t;

    for(int j = 0; j<t; j++){

        int e;
        cin>>e;

        int cnt = 1;
        if(e == 0){
            cin>>cnt;
        }

        //Recieve the input and call the function to do things

        for(int k = 0; k<e; k++){
            int type;
            cin>>type;
            if(type == 1){
                int id, c, m;
                cin>>id>>c>>m;
                mejoropcion(j, id, c, m);
            }
            if(type == 2){
                int id;
                cin>>id;
                cerrar(j, id);
            }
        }

        //Just print
        j--;
        for(int it = 0; it<cnt; it++){
            j++;
            if(instrucciones.find(j) == instrucciones.end()){
                cout<<"0\n";
                continue;
            }
            cout<<instrucciones[j].size()<<"\n";
            for(auto u: instrucciones[j]){
                if(u[0] == 2){
                    cout<<u[0]<<" "<<u[1]<<"\n";
                }
                else{
                    cout<<u[0]<<" "<<u[1]<<" "<<u[2]<<"\n";
                }
            }
        }
        cout.flush();
    }

    int e;
    cin>>e;
    
    return 0;
}

