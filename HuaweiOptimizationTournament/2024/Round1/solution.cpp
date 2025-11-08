#pragma GCC optimize("Ofast,unroll-loops")
#include<bits/stdc++.h>
#if defined(_WIN32) || defined(_WIN64)
#include <Windows.h>
#endif
using namespace std;
typedef long long ll;
typedef pair<int, int> ii;

//Define the global variables
const ll BILLION = 1e9;
ll l, h, F;
double Fd, ld, hd;
double alp, bet, gam;
int n, m;
vector<vector<double> > fI;
vector<vector<double> > fO;
vector<vector<double> > fIO;
vector<vector<int> > g;
vector<int> Pres;
vector<int> N;
vector<double> tau;
vector<int> candidato = {0, 0, 0};
double Lopt, Lpopt, Ldopt;
double Lmax, Lpmax, Ldmax;
double umax = 0;
double fmaxx = 0;
double Nmin = 1e18;
double IOmin = 1e18;
double Imin = 1e18;
double Omin = 1e18;
double lh8;
ll checkA, checkB;
mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());
 
//Function to get the time 
long double cpu_time() {
#if defined(_WIN32) || defined(_WIN64)
    FILETIME creation_ft, exit_ft, kernel_ft, user_ft;
    GetProcessTimes(GetCurrentProcess(), &creation_ft, &exit_ft, &kernel_ft, &user_ft);
 
    auto extract_time = [](FILETIME ft) {
        return 1e-7L * (ft.dwLowDateTime | uint64_t(ft.dwHighDateTime) << 32);
    };
 
    return extract_time(user_ft) + extract_time(kernel_ft);
#else
    return (long double) clock() / CLOCKS_PER_SEC;
#endif
}

void recalculate();

//Class that represent a pipeline
struct pipeline{

    //Define the variables for each pipeline
    bool activo = true;
    int index;
    double t, f, c, e, sumL, sumLp, sumLd, ft;
    double savesumL, savesumLp, savesumLd;
    vector<double> L, Lp, Ldmc, Lc, v;
    int jprimero = -1;
    int rprimero = -1;
    
    //The init, we start the variables
    void init(ll _t, ll _f, ll _c, ll _e, int _index){
        index = _index;
        t = _t;
        f = _f;
        c = _c;
        e = _e;
        ft = f*t;
        sumL = 0;
        for(int j = 0; j<m; j++){
            sumL += tau[j];
        }
        sumLp = 0;
        sumLd = 0;
        L.resize(m);
        Lp.resize(m);
        Ldmc.resize(m);
        Lc.resize(m);
        v.resize(m);
    }

    //We add the request (j, r) to this pipeline
    void agregar(int j, int r){

        g[j][r] = index;

        sumL -= max(Lp[j] + Ldmc[j] + Lc[j], tau[j]);
        sumLd -= Ldmc[j];

        Lp[j] += (fI[j][r] + fO[j][r])/ft;
        sumLp += fI[j][r]/ft;
        sumLd += fO[j][r]/ft;
        v[j] += fIO[j][r];
        Ldmc[j] = (Fd)/(t*c) + (lh8*v[j])/(t*c);
        Lc[j] = (lh8*v[j]*(t - 1.0))/(e*t);

        sumL += max(Lp[j] + Ldmc[j] + Lc[j], tau[j]);
        sumLd += Ldmc[j];

        Lmax = max(Lmax, sumL);
        Lpmax = max(Lpmax, sumLp);
        Ldmax = max(Ldmax, sumLd);

        if(jprimero == -1){
            jprimero = j;
            rprimero = r;
        }

        if(sumL == Lmax){
            candidato[0] = index;
        }
        if(sumLp == Lpmax){
            candidato[1] = index;
        }
        if(sumLd == Ldmax){
            candidato[2] = index;
        }
    }

    //we add virtually the request (j, r) to this pipeline
    void xagregar(int j, int r){

        sumL -= max(Lp[j] + Ldmc[j] + Lc[j], tau[j]);
        sumLd -= Ldmc[j];

        Lp[j] += (fI[j][r] + fO[j][r])/ft;
        sumLp += fI[j][r]/ft;
        sumLd += fO[j][r]/ft;
        v[j] += fIO[j][r];
        Ldmc[j] = (Fd)/(t*c) + (lh8*v[j])/(t*c);
        Lc[j] = (lh8*v[j]*(t - 1.0))/(e*t);

        sumL += max(Lp[j] + Ldmc[j] + Lc[j], tau[j]);
        sumLd += Ldmc[j];
    }

    //we erase the request (j, r) to this pipeline
    void quitar(int j, int r){
        
        bool recalcular = ((sumL == Lmax) || (sumLp == Lpmax) || (sumLd == Ldmax));

        g[j][r] = -1;
        sumL -= max(Lp[j] + Ldmc[j] + Lc[j], tau[j]);
        sumLd -= Ldmc[j];

        Lp[j] -= (fI[j][r] + fO[j][r])/(f*t);
        sumLp -= fI[j][r]/(f*t);
        sumLd -= fO[j][r]/(f*t);
        v[j] -= fIO[j][r];
        Ldmc[j] = (Fd)/(t*c) + (lh8*v[j])/(t*c);
        Lc[j] = (lh8*v[j]*(t - 1.0))/(e*t);

        sumL += max(Lp[j] + Ldmc[j] + Lc[j], tau[j]);
        sumLd += Ldmc[j];

        if(recalcular == true){
            recalculate();
        }

        if(j == jprimero){
            jprimero = -1;
        }

    }

    //we erase virtually the request (j, r) to this pipeline
    void xquitar(int j, int r){

        sumL -= max(Lp[j] + Ldmc[j] + Lc[j], tau[j]);
        sumLd -= Ldmc[j];

        Lp[j] -= (fI[j][r] + fO[j][r])/(f*t);
        sumLp -= fI[j][r]/(f*t);
        sumLd -= fO[j][r]/(f*t);
        v[j] -= fIO[j][r];
        Ldmc[j] = (Fd)/(t*c) + (lh8*v[j])/(t*c);
        Lc[j] = (lh8*v[j]*(t - 1.0))/(e*t);

        sumL += max(Lp[j] + Ldmc[j] + Lc[j], tau[j]);
        sumLd += Ldmc[j];
    }

};

//We initialize the pipelines
vector<pipeline> P;

//Recalculate Lmax, Lpmax and Ldmax to get the score function
void recalculate(){
    Lmax = 0;
    Lpmax = 0;
    Ldmax = 0;
    for(int i = 0; i<int(P.size()); i++){
        Lmax = max(Lmax, P[i].sumL);
        Lpmax = max(Lpmax, P[i].sumLp);
        Ldmax = max(Ldmax, P[i].sumLd);
    }
}

//Return a relaxation of the score function
double getscore(){
    return alp*Lopt/Lmax + bet*Lpopt/Lpmax + gam*Ldopt/Ldmax;
}

//Return the score function
double realscore(){
    return floor(alp*(floor(1e7*Lopt/Lmax)) + bet*(floor(1e7*Lpopt/Lpmax)) + gam*(floor(1e7*Ldopt/Ldmax)));
}

//Output the answer in the format said by the problem
void print(){

    for(int i = 0; i<n; i++){
        cout<<Pres[i]<<" "<<8/Pres[i]<<" 1\n";
    }

    int pcount = 0;
    vector<int> realindex(P.size());
    for(int i = 0; i<int(P.size()); i++){
        if(P[i].activo == false) continue;
        realindex[i] = pcount;
        pcount++;
    }

    for(int j = 0; j<m; j++){
        vector<int> W(P.size(), 0);
        for(int r = 0; r<N[j]; r++){
            assert(P[g[j][r]].activo == true);
            cout<<realindex[g[j][r]] + 1<<" "<<W[g[j][r]] + 1<<"\n";
            W[g[j][r]]++;
        }
    }
}

//Try a random change of requests (j1, r1) and (j2, r2) between respective pipelines
void cambiar(bool samej, double temperatura){

    int j1 = rng()%m;
    int r1 = rng()%N[j1];

    int j2 = j1;
    if(samej == false){
        j2 = rng()%m;
    }
    
    int r2 = rng()%N[j2];
    
    int i1 = g[j1][r1];
    int i2 = g[j2][r2];

    if(i1 == i2) return;

    double antes = getscore();

    P[i1].quitar(j1, r1);
    P[i2].quitar(j2, r2);

    P[i1].agregar(j2, r2);
    P[i2].agregar(j1, r1);

    double delta = (getscore() - antes)/antes;

    double random_number = static_cast<double>(rand())/RAND_MAX;

    if(delta >= 0 || random_number <= exp(delta/temperatura)){
        return;
    }

    P[i1].quitar(j2, r2);
    P[i2].quitar(j1, r1);

    P[i1].agregar(j1, r1);
    P[i2].agregar(j2, r2);
}

//Put out requests from the 3 pipelines that are actives as maximum of some latency
//and assign it to a new pipeline greedyly
void sacar(){
    queue<ii> q;
    if(P[candidato[0]].jprimero != -1){
        q.push({P[candidato[0]].jprimero, P[candidato[0]].rprimero});
    }
    if(P[candidato[1]].jprimero != -1){
        q.push({P[candidato[1]].jprimero, P[candidato[1]].rprimero});
    }
    if(P[candidato[2]].jprimero != -1){
        q.push({P[candidato[2]].jprimero, P[candidato[2]].rprimero});
    }
    while(q.size() > 0){
        int j = q.front().first; int r = q.front().second;
        q.pop();
        P[g[j][r]].quitar(j, r);
        double maxi = 0;
        int where = -1;
        int shift = rng()%P.size();
        for(int ix = 0; ix<int(P.size()); ix++){
            int i = (ix + shift)%P.size();
            if(P[i].activo == false) continue;
            P[i].xagregar(j, r);
            double score = alp*Lopt/max(Lmax, P[i].sumL) + bet*Lpopt/max(Lpmax, P[i].sumLp) + gam*Ldopt/max(Ldmax, P[i].sumLd);
            if(score >= maxi){
                maxi = score;
                where = i;
            }
            P[i].xquitar(j, r);
        }
        P[where].agregar(j, r);
    }
}

//We try to not put the maximum possible pipelines
//And if the answer is better, we do it
bool achicar(int unit, int pstart){
    if(Pres[unit] == 1) return false;
    double scoreantes = getscore();
    queue<ii> q;
    for(int j = 0; j<m; j++){
        for(int r = 0; r<N[j]; r++){
            if(g[j][r] >= pstart && g[j][r] < pstart + Pres[unit]){
                P[g[j][r]].xquitar(j, r);
                q.push({j, r});
            }
        }
    }
    for(int i = pstart; i<pstart + Pres[unit]/2; i++){
        P[i].t = 2.0*P[i].t;
        P[i].ft = 2.0*P[i].ft;
    }
    recalculate();
    while(q.size() > 0){
        int j = q.front().first;
        int r = q.front().second;
        q.pop();
        double maxi = 0;
        int where = -1;
        for(int i = pstart; i<pstart + Pres[unit]/2; i++){
            P[i].xagregar(j, r);
            double score = alp*Lopt/max(Lmax, P[i].sumL) + bet*Lpopt/max(Lpmax, P[i].sumLp) + gam*Ldopt/max(Ldmax, P[i].sumLd);
            if(score >= maxi){
                maxi = score;
                where = i;
            }
            P[i].xquitar(j, r);
        }
        P[where].agregar(j, r);
    }
    
    if(getscore() >= scoreantes){
        for(int i = pstart + Pres[unit]/2; i<pstart+Pres[unit]; i++){
            P[i].activo = false;
        }
        Pres[unit] = Pres[unit]/2;
        return true;
    }
    else{
        for(int j = 0; j<m; j++){
            for(int r = 0; r<N[j]; r++){
                if(g[j][r] >= pstart && g[j][r] < pstart + Pres[unit]){
                    P[g[j][r]].xquitar(j, r);
                    q.push({j, r});
                }
            }
        }
        recalculate();
        for(int i = pstart; i<pstart + Pres[unit]; i++){
            P[i].t = 8.0/double(Pres[unit]);
            P[i].ft = P[i].t*P[i].f;
        }
        while(q.size() > 0){
            int j = q.front().first;
            int r = q.front().second;
            q.pop();
            double maxi = 0;
            int where = -1;
            for(int i = pstart; i<pstart + Pres[unit]; i++){
                P[i].xagregar(j, r);
                double score = alp*Lopt/max(Lmax, P[i].sumL) + bet*Lpopt/max(Lpmax, P[i].sumLp) + gam*Ldopt/max(Ldmax, P[i].sumLd);
                if(score >= maxi){
                    maxi = score;
                    where = i;
                }
                P[i].xquitar(j, r);
            }
            P[where].agregar(j, r);
            Lmax = max(Lmax, P[where].sumL);
            Lpmax = max(Lpmax, P[where].sumLp);
            Ldmax = max(Ldmax, P[where].sumLd);
        }
        return false;
    }
}

int main(){
	
    ios::sync_with_stdio(false);
    cin.tie(0);

    //Read the input
	cin>>l>>h>>F;
    Fd = 2LL*F; ld = l; hd = h;
	cin>>alp>>bet>>gam;
	cin>>n>>m;

    vector<vector<ll> > units(n, vector<ll>(5));

	for(int i = 0; i<n; i++){
        cin>>units[i][0]>>units[i][1]>>units[i][2]>>units[i][3]>>units[i][4];
        units[i][1] *= BILLION;
		units[i][2] *= BILLION;
		units[i][3] *= BILLION;
		units[i][4] *= BILLION;
        umax = max(umax, double(units[i][0]));
        fmaxx = max(fmaxx, double(units[i][1]));
	}

    //Resize the global variables and assign the values to it
    N.resize(m);
    tau.resize(m, 0);
    fI.resize(m);
    fO.resize(m);
    fIO.resize(m);
    g.resize(m);
    vector<double> realtau(m);
    ll IOmax = 0;

    for(int j = 0; j<m; j++){
        cin>>N[j]>>realtau[j];
        fI[j].resize(N[j]);
        fO[j].resize(N[j]);
        fIO[j].resize(N[j]);
        g[j].resize(N[j]);
        Nmin = min(Nmin, double(N[j]));
        vector<double> I(N[j]);
        for(int r = 0; r<N[j]; r++){
            cin>>I[r];
            Imin = min(Imin, I[r]);
            fI[j][r] = Fd*I[r];
        }
        for(int r = 0; r<N[j]; r++){
            double O;
            cin>>O;
            fO[j][r] = Fd*O;
            fIO[j][r] = O*(I[r] + (O - 1.0)/2.0);
            Omin = min(Omin, O);
            IOmin = min(IOmin, I[r] + O);
            IOmax = max(IOmax, ll(I[r]) + ll(O));
        }
    }

    checkA = 2LL*F;
    checkB = 4LL*l*h*IOmax;
    Pres.resize(n);
    //Create as many pipelines as possible
    for(int i = 0; i<n; i++){
        ll p = 1; ll t = 8;
        for(auto _p: {1, 2, 4, 8}){
            ll _t = units[i][0]/_p;
            if(units[i][2]*_t >= checkA + checkB){
                p = _p;
                t = _t;
            }
        }
        Pres[i] = p;
        for(int j = 0; j<p; j++){
            pipeline pipe;
            pipe.init(t, units[i][1], units[i][3], units[i][4], P.size());
            P.push_back(pipe);
        }
    }

    //Calculate constant to get the score in the future
    Lopt = (Fd*double(m)*Nmin*IOmin)/(umax*umax*fmaxx);
    Lpopt = (Fd*double(m)*Nmin*Imin)/(umax*umax*fmaxx);
    Ldopt = (Fd*double(m)*Nmin*Omin)/(umax*umax*fmaxx);
    lh8 = 8.0*ld*hd;

    Lmax = 0;
    Lpmax = 0;
    Ldmax = 0;

    //We put each request (j, r) greedily
    //And assuming tau = 0
    for(int j = 0; j<m; j++){
        for(int r = 0; r<N[j]; r++){
            double maxi = 0;
            int where = -1;
            int shift = rng()%P.size();
            for(int ix = 0; ix<int(P.size()); ix++){
                int i = (ix + shift)%P.size();
                P[i].xagregar(j, r);
                double score = alp*Lopt/max(Lmax, P[i].sumL) + bet*Lpopt/max(Lpmax, P[i].sumLp) + gam*Ldopt/max(Ldmax, P[i].sumLd);
                if(score >= maxi){
                    maxi = score;
                    where = i;
                }
                P[i].xquitar(j, r);
            }
            g[j][r] = where;
            P[where].agregar(j, r);
            Lmax = max(Lmax, P[where].sumL);
            Lpmax = max(Lpmax, P[where].sumLp);
            Ldmax = max(Ldmax, P[where].sumLd);
        }
    }

    //We assign the real value for tau
    for(int j = 0; j<m; j++){
        tau[j] = realtau[j];
    }

    //Recalculate the L for each pipeline with the new tau
    for(int i = 0; i<int(P.size()); i++){
        for(int j = 0; j<m; j++){
            P[i].sumL -= (P[i].Lp[j] + P[i].Ldmc[j] + P[i].Lc[j]);
            P[i].sumL += max(P[i].Lp[j] + P[i].Ldmc[j] + P[i].Lc[j], tau[j]);
        }
    }

    //Recalculate the max latency
    recalculate();

    //We try to not put as many pipelines as possible
    int pstart = 0;
    vector<int> Presx = Pres;
    for(int i = 0; i<n; i++){
        while(achicar(i, pstart) == true);
        pstart += Presx[i];
    }

    //Get the time
    double timestart = cpu_time();
    double timefinish = 3.9;
    double delta = timefinish - timestart;

    //Do annealing while we have time remaining
    while(cpu_time() < timefinish){
        double temperatura = exp(-9.0 - 18.0*(cpu_time() - timestart)/delta);
        //We try 45 local changes
        for(int colocolo = 0; colocolo < 45; colocolo++){
            cambiar(true, temperatura);
            cambiar(false, temperatura);
        }
        //We try 1 'global' change
        sacar();
    }

    //Print the answer
    print();

    return 0;
}