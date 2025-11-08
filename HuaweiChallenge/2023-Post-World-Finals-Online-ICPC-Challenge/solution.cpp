#pragma GCC optimize("-fexcess-precision=fast")
#pragma GCC target("fpmath=387")
#include<bits/stdc++.h>
#if defined(_WIN32) || defined(_WIN64)
#include <Windows.h>
#endif
using namespace std;
const double INF = std::numeric_limits<double>::infinity();
typedef pair<int, int> ii;
typedef long long ll;
typedef unsigned long long ull;
const ull mascara = (((~(ull(1)<<ull(63)))<<ull(53))>>ull(1));
const ull mascara2 = ((~(ull(0)))<<ull(42));
const ull mascara3 = (ull(1039)<<ull(52));
const ull mascara4 = (ull(1008)<<ull(52));
const vector<vector<int> > recpenalty = {{}, {0, 0, -2, -3}, {0, 2, 0, -1}, {0, 3, 1, 0}};
const vector<vector<int> > opcionesp = {{}, {-2, -2}, {-1, -3}, {-2, -2}};
mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());

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

inline double fastPow(double a, double b){
    union{
        double d;
        int x[2];
    }u = {a};
    u.x[1] = (int)(b * (u.x[1] - 1072632447) + 1072632447);
    u.x[0] = 0;
    return u.d;
}

double realsum(const vector<double>& vec){
    long double trueSum=0, corr=0;
    vector<double> dvtmp=vec;
    sort(dvtmp.begin(),dvtmp.end(), [](const double x, const double y) {
        return fabs(x) < fabs(y);
    });
    for (auto i : dvtmp) {
        volatile long double y = static_cast<long double>(i) - corr;
        volatile long double t = trueSum + y;
        corr = (t - trueSum) - y;
        trueSum = t;
    }
    return (double)trueSum;
}

double sumd(double a, double b){
    return a + b;
}

double sums(double a, double b){
    a = float(a);
    b = float(b);
    return a + b;
}

double tohalf(double a){
    ull bits = *reinterpret_cast<ull*>(&a);
    if((bits&mascara) > mascara3){
        return INF;
    }
    if((bits&mascara) < mascara4){
        return 0.0;
    }   
    bits = bits&mascara2;
    a = *reinterpret_cast<double*>(&bits);
    return a;
}

double sumh(double a, double b){
    a = tohalf(a);
    b = tohalf(b);
    return tohalf(a + b);
}

struct Arbol{

    int n, N, penalty, penaltysave;
    double real, fake, score, bestscore, distancia, den1, den2;
    bool parchado = false;
    vector<double> v, vsave;
    vector<int> padre, cambioL, cambioR, cambioP, padresave;
    vector<ii> hijos, hijossave;
    vector<int> type, typesave;

    Arbol(vector<pair<double, int> > &a, double reall){

        n = a.size();
        N = 2*n - 1;
        penalty = 4*(n-1);
        real = reall;
        type.resize(N, -1);
        v.resize(N);
        hijos.resize(N, {-1, -1});
        padre.resize(N, -1);
        cambioL.resize(N);
        cambioR.resize(N);
        vector<ii> borde(N);
        bestscore = -10;
        den1 = double(n-1);
        den2 = max(abs(real), 1e-200);

        queue<int> q;
        borde[0].first = 0;
        borde[0].second = n-1;
        int ahora = 1;
        q.push(0);
        
        while(q.size() > 0){
            int nodo = q.front();
            q.pop();
            if(borde[nodo].first == borde[nodo].second){
                v[nodo] = a[borde[nodo].first].first;
                type[nodo] = a[borde[nodo].first].second;
                continue;
            }
            cambioP.push_back(nodo);
            int largo = borde[nodo].second - borde[nodo].first + 1;
            int largol;
            if(largo > 16){
                largol = 16*((largo + 15)/32);
            }
            else{
                largol = largo/2;
            }
            hijos[nodo].first = ahora;
            borde[ahora].first = borde[nodo].first;
            borde[ahora].second = borde[nodo].first + largol - 1;
            padre[ahora] = nodo;
            q.push(ahora);
            ahora++;
            hijos[nodo].second = ahora;
            borde[ahora].first = borde[nodo].first + largol;
            borde[ahora].second = borde[nodo].second;
            padre[ahora] = nodo;
            q.push(ahora);
            ahora++;
        }

        bool cortado = false;
        if(n%16 != 0){
            cortado = true;
        }

        cambioL[0] = 0;

        for(int i = 1; i<N; i++){
            cambioL[i] = i;
            int largo = borde[i].second - borde[i].first;
            if(largo <= 16){
                if(borde[i].first == borde[i-1].second + 1 && borde[i-1].first/16 == borde[i].first/16){
                    cambioL[i] = cambioL[i-1];
                }
                continue;
            }
            if(borde[i].second == n-1 && cortado == true){
                continue;
            }
            if(borde[i].first == borde[i-1].second + 1){
                cambioL[i] = cambioL[i-1];
            }
        }

        cambioR[N-1] = N-1;

        for(int i = N-2; i>=0; i--){
            cambioR[i] = i;
            int largo = borde[i].second - borde[i].first;
            if(largo <= 16){
                if(borde[i].second + 1 == borde[i+1].first && borde[i].first/16 == borde[i+1].first/16){
                    cambioR[i] = cambioR[i+1];
                }
                continue;
            }
            if(borde[i+1].second == n-1 && cortado == true){
                continue;
            }
            if(borde[i].second + 1 == borde[i+1].first){
                cambioR[i] = cambioR[i+1];
            }
        }

        for(int i = N-1; i>=0; i--){
            if(type[i] < 0){
                sumar(i);
            }
        }
        save();
        distancia = 1;
        computescore();
    }

    void change(int i, int j){
        
        if(hijos[padre[i]].first == i){
            hijos[padre[i]].first = j;
        }
        else{
            hijos[padre[i]].second = j;
        }

        if(hijos[padre[j]].first == j){
            hijos[padre[j]].first = i;
        }
        else{
            hijos[padre[j]].second = i;
        }

        swap(padre[i], padre[j]);

        do{
            if(i > j){
                i = padre[i];
                if(sumar(i)) i = 0;
            }
            if(j > i){
                j = padre[j];
                if(sumar(j)) j = 0;
            }
        }while(i != j);

        while(i != -1){
            sumar(i);
            i = padre[i];
        }

        computescore();
    }

    void ejecutarchange(){
        int i = rng()%N;
        i = rng()%(i+1);
        i = rng()%(i+1);
        i = rng()%(i+1);
        int j = rng()%(cambioR[i] - cambioL[i] + 1) + cambioL[i];
        if(i == j) return;
        double scoreantes = score;
        change(i, j);
        if(scoreantes > score){
            change(i, j);
        }
    }

    void precision(int i, int nuevo){
        penalty = penalty + recpenalty[-type[i]][-nuevo];
        type[i] = nuevo;
        while(i != -1){
            if(sumar(i)) break;
            i = padre[i];
        }
        computescore();
    }

    void ejecutarprecision(){
        double scoreantes = score;
        int i = rng()%(cambioP.size());
        int tipoantes = type[cambioP[i]];
        int nuevo = opcionesp[-type[cambioP[i]]][rng()%2];
        precision(cambioP[i], nuevo);
        if(scoreantes > score){
            precision(cambioP[i], tipoantes);
        }
    }

    bool precisionb(int i, int nuevo){
        double scoreantes = score;
        penalty = penalty + recpenalty[-type[i]][-nuevo];
        type[i] = nuevo;
        while(i != -1){
            double vantes = v[i];
            if(sumar(i)) break;
            if((fake < real && v[i] < vantes) || (fake > real && v[i] > vantes)){
                return false;
            }
            i = padre[i];
        }
        computescore();
        if(scoreantes > score){
            return false;
        }
        return true;
    }

    void ejecutarprecisionb(){
        int i = rng()%(cambioP.size());
        int tipoantes = type[cambioP[i]];
        int nuevo = opcionesp[-type[cambioP[i]]][rng()%2];
        if(precisionb(cambioP[i], nuevo) == false){
            precision(cambioP[i], tipoantes);
        }
    }

    void posinicial(double cota1, double cota2, int x = 0, double multiplicador = 1.45){
        if(type[x] >= 0) return;
        posinicial(multiplicador*cota1, multiplicador*cota2, hijos[x].first);
        posinicial(multiplicador*cota1, multiplicador*cota2, hijos[x].second);
        double d = v[hijos[x].first] + v[hijos[x].second];
        double s = float(v[hijos[x].first]) + float(v[hijos[x].second]);
        double h = tohalf(tohalf(v[hijos[x].first]) + tohalf(v[hijos[x].second]));
        int termino = type[x];
        if(isfinite(s) && abs(s - d) <= cota1){
            termino = min(termino, -2);
        }
        if(isfinite(h) && abs(h - d) <= cota2){
            termino = min(termino, -3);
        }
        penalty = penalty + recpenalty[-type[x]][-termino];
        type[x] = termino;
        //precision(x, termino);
        sumar(x);
    }

    void binary(){
        double inf1 = 1e-30;
        double sup1 = 1;
        while(cpu_time() < 3.5){
            double mid = sqrt(inf1*sup1);
            posinicial(1e-30, mid);
            computescore();
            for(int j = 0; j<400000; j++){
                ejecutarchange();
                ejecutarprecisionb();
                if(distancia <= 1e-20) break;
            }
            if(score >= bestscore || distancia <= 1e-20){
                save();
                bestscore = score;
                inf1 = mid;
            }
            else{
                load();
                sup1 = mid;
            }
        }
        load();
        double inf2 = 1e-30;
        double sup2 = inf1;
        while(cpu_time() < 6.0){
            double mid = sqrt(inf2*sup2);
            posinicial(mid, inf1);
            computescore();
            for(int j = 0; j<400000; j++){
                ejecutarchange();
                ejecutarprecisionb();
                if(distancia <= 1e-20) break;
            }
            if(score >= bestscore || distancia <= 1e-20){
                save();
                bestscore = score;
                inf2 = mid;
            }
            else{
                load();
                sup2 = mid;
            }
        }
        load();
        if(score < 20 && n <= 200000){
            parche();
        }
        if(parchado == true){
            return;
        }
        while(cpu_time() < 9.5){
            for(int i = 0; i<1000; i++){
                ejecutarchange();
                ejecutarprecision();
            }
        }
    }

    void parche(){
        penalty = 0;
        set<pair<double, int> > S;
        queue<int> q;
        for(int i = N-1; i>=0; i--){
            if(type[i] >= 0){
                S.insert({v[i], i});
            }
            else{
                q.push(i);
            }
        }
        while(S.size() > 1){
            auto it1 = S.begin();
            pair<double, int> uno = *it1;
            S.erase(it1);
            auto it2 = S.end();
            it2--;
            pair<double, int> dos = *it2;
            S.erase(it2);
            int nodo = q.front();
            q.pop();
            hijos[nodo].first = uno.second;
            hijos[nodo].second = dos.second;
            double d = uno.first + dos.first;
            double s = float(uno.first) + float(dos.first);
            double h = tohalf(tohalf(uno.first) + tohalf(dos.first));
            v[nodo] = d;
            if(h == d){
                type[nodo] = -3;
                penalty = penalty + 1;
                S.insert({d, nodo});
                continue;
            }
            if(s == d){
                type[nodo] = -2;
                penalty = penalty + 2;
                S.insert({d, nodo});
                continue;
            }
            type[nodo] = -1;
            penalty = penalty + 4;
            S.insert({d, nodo});
        }
        double scoreantes = score;
        computescore();

        if(score > scoreantes){
            parchado = true;
        }
        else{
            load();
        }
    }

    void recomputepenalty(int i, int nuevo){
        penalty = penalty + recpenalty[-type[i]][-nuevo];
        type[i] = nuevo;
    }

    bool sumar(int i){
        if(type[i] == -1){
            return v[i] == (v[i] = v[hijos[i].first] + v[hijos[i].second]);
        }
        if(type[i] == -2){
            return v[i] == (v[i] = float(v[hijos[i].first]) + float(v[hijos[i].second]));
        }
        return v[i] == (v[i] = tohalf(tohalf(v[hijos[i].first]) + tohalf(v[hijos[i].second])));
    }

    void load(){
        penalty = penaltysave;
        hijos = hijossave;
        padre = padresave;
        v = vsave;
        type = typesave;
        computescore();
    }

    void save(){
        penaltysave = penalty;
        hijossave = hijos;
        padresave = padre;
        vsave = v;
        typesave = type;
    }
    
    void computescore(){
        fake = v[0];
        if(isfinite(fake) == false){
            score = -10;
            distancia = 1e100;
            return;
        }
        distancia = max(1e-20, abs(real - fake)/den2);
        double D = 10.0/sqrt(0.5 + double(penalty)/den1);
        double A = fastPow(distancia, 0.05);
        score = D/A;
    }

    void imprimir(string &s, int x = 0){
        if(type[x] >= 0){
            for(auto ch: to_string(type[x] + 1)){
                s.push_back(ch);
            }
            return;
        }
        s.push_back('{');
        if(type[x] == -1){
            s.push_back('d');
        }
        if(type[x] == -2){
            s.push_back('s');
        }
        if(type[x] == -3){
            s.push_back('h');
        }
        s.push_back(':');
        imprimir(s, hijos[x].first);
        s.push_back(',');
        imprimir(s, hijos[x].second);
        s.push_back('}');
        return;
    }
    
};

int main(){

    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n;
    cin>>n;

    vector<pair<double, int> > A(n);
    vector<double> a(n);

    for(int i = 0; i<n; i++){
        string x;
        cin>>x;
        a[i] = stod(x);
    }

    for(int i = 0; i<n; i++){
        A[i].first = a[i];
        A[i].second = i;
    }

    double real = realsum(a);
    a.clear();

    vector<pair<double, int> > grande;
    for(int i = 0; i<n; i = i + 16){
        sort(A.begin() + i, A.begin() + min(i+16, n), [](const pair<double, int> x, const pair<double, int> y) {return abs(x.first) < abs(y.first);});
        double sum = 0;
        int sum2 = 0;
        for(int j = 0; j<16 && i+j<n; j++){
            sum = sum + A[i+j].first;
            sum2++;
        }
        if(sum2 == 16){
            grande.push_back({sum, i/16});
        }
    }
    sort(grande.begin(), grande.begin(), [](const pair<double, int> x, const pair<double, int> y) {return abs(x.first) < abs(y.first);});
    
    vector<pair<double, int> > Ax = A;

    for(int i = 0; i<int(grande.size()); i++){
        for(int k = 0; k<16; k++){
            Ax[16*i + k] = A[16*grande[i].second + k];
        }
    }

    Arbol arbol(Ax, real);

    arbol.binary();
    string res = "";
    arbol.imprimir(res);
    cout<<res<<endl;
    cerr<<arbol.score<<endl;

    return 0;
}//f10
