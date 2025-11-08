#pragma GCC optimize("Ofast")
#include<bits/stdc++.h>
using namespace std;
typedef pair<int, int> ii;
int N, K, T, R, J;
const double W = 192/log(2);
const double duno = 1;
vector<vector<vector<vector<double> > > > s0;
vector<vector<vector<vector<double> > > > d;
vector<vector<vector<vector<double> > > > p;
vector<vector<int> > A;

double parseador(string &s){
    double res = 0;
    double signo = 1;
    double aux = 0.1;
    bool cambio = false;
    for(auto ch: s){
        if(ch == '-'){
            signo = -1;
            continue;
        }
        if(ch == '.'){
            cambio = true;
            continue;
        }
        if(cambio == false){
            res = 10.0*res + double(ch - '0');
        }
        else{
            res = res + aux*(double(ch - '0'));
            aux = aux*0.1;
        }
    }
    return signo*res;
}

int LocalMaxSinInterUpdate(int indicesort){

    for(int t = 0; t<T; t++){
        for(int k = 0; k<K; k++){
            for(int r = 0; r<R; r++){
                for(int i = 0; i<N; i++){
                    p[t][i][k][r] = 0;
                }
            }
        }
    }

    int cuantos = 0;
    for(int i = 0; i<J; i++){
        A[i][0] = A[i][indicesort];
    }

    sort(A.begin(), A.end());

    vector<vector<vector<double> > > gananciax(J, vector<vector<double> >(T, vector<double>(R, 0)));
    vector<set<pair<double, vector<int> > > > ganancia(J);
    vector<double> libre(J, 0);
    vector<double> suma(J, 0);
    vector<double> cota(J);
    vector<vector<int> > quien(T, vector<int>(R, -1));

    for(int j = 0; j<J; j++){
        cota[j] = A[j][2];
        for(int t = A[j][4]; t<=A[j][4] + A[j][5] - 1; t++){
            for(int r = 0; r<R; r++){
                double mean = 1;
                for(int k = 0; k<K; k++){
                    mean = mean*pow(s0[t][A[j][3]][k][r], duno/double(K));
                }
                gananciax[j][t][r] = W*K*log1p(mean);
                ganancia[j].insert({-W*K*log1p(mean), {t, r}});
                libre[j] = libre[j] + W*K*log1p(mean);
            }
        }
    }

    for(int j = 0; j<J; j++){
        int id = A[j][3];
        
        while(ganancia[j].size() > 0){
            if(suma[j] >= cota[j]) break;
            
            auto nodo = *ganancia[j].begin();
            ganancia[j].erase(nodo);
            double profit = -nodo.first;
            int t = nodo.second[0];
            int r = nodo.second[1];
            int exid = quien[t][r];
            
            if(exid == -1){
                suma[j] = suma[j] + profit;
                quien[t][r] = j;
                for(int k = 0; k<K; k++){
                    p[t][id][k][r] = 1;
                }
                for(int jx = 0; jx<J; jx++){
                    libre[jx] = libre[jx] - gananciax[jx][t][r];
                }
                continue;
            }
            
            if(exid != -1 && libre[exid] + (suma[exid] - cota[exid]) - gananciax[exid][t][r] >= 0){
                suma[j] = suma[j] + profit;
                quien[t][r] = j;
                for(int k = 0; k<K; k++){
                    p[t][id][k][r] = 1;
                    p[t][A[exid][3]][k][r] = 0;
                }
                suma[exid] = suma[exid] - gananciax[exid][t][r];
                while(ganancia[exid].size() > 0){
                    if(suma[exid] >= cota[exid]) break;
                    auto nodo2 = *ganancia[exid].begin();
                    ganancia[exid].erase(nodo2);
                    if(quien[nodo2.second[0]][nodo2.second[1]] != -1) continue;
                    for(int k = 0; k<K; k++){
                        p[nodo2.second[0]][A[exid][3]][k][nodo2.second[1]] = 1;
                    }
                    suma[exid] = suma[exid] - nodo2.first;
                    quien[nodo2.second[0]][nodo2.second[1]] = exid;
                    for(int jx = 0; jx<J; jx++){
                        libre[jx] = libre[jx] - gananciax[jx][nodo2.second[0]][nodo2.second[1]];
                    }
                }
            }
        }
        
        if(suma[j] < cota[j]){
            for(int t = A[j][4]; t<=A[j][4] + A[j][5] - 1; t++){
                for(int r = 0; r<R; r++){
                    if(quien[t][r] == j){
                        quien[t][r] = -1;
                        for(int k = 0; k<K; k++){
                            p[t][id][k][r] = 0;
                        }
                        for(int jx = 0; jx<J; jx++){
                            if(t >= A[jx][4] && t <= A[jx][4] + A[jx][5] - 1){
                                libre[jx] = libre[jx] + gananciax[jx][t][r];
                                ganancia[jx].insert({-gananciax[jx][t][r], {t, r}});
                            }
                        }
                    }
                }
            }
        }
        else{
            cuantos = cuantos + 1;
        }
    }

    return cuantos;
}

pair<int, vector<vector<vector<double> > > > Tunosolve(vector<pair<double, int> > &Q, int t, int iteraciones = 500){
//O(iteraciones*n^2*r)
    sort(Q.begin(), Q.end());

    int n = min(int(Q.size()), 50);

    vector<int> bitmask(R, 0);
    int maxi = 0;
    vector<int> bitmaskmax(R, 0);
    vector<vector<double> > pauxmaxi;
    vector<vector<double> > paux(n, vector<double>(R, 0));
    vector<double> sobra(R);
    vector<bool> fuera(n);
    vector<int> cuantos(K);
    vector<double> cota(n);
    double dRx = double(min(4, R));

    for(int xd = 0; xd<iteraciones; xd++){

        int sum = 0;

        for(int k = 0; k<K; k++){
            cuantos[k] = 0;
        }

        for(int r = 0; r<R; r++){
            bitmask[r] = rand()%K;
            sobra[r] = 0;
            cuantos[bitmask[r]] += n;
        }

        for(int i = 0; i<n; i++){
            fuera[i] = false;
            cota[i] = 0;
            for(int r = 0; r<R; r++){
                paux[i][r] = dRx/double(cuantos[bitmask[r]]);
            }
        }

        for(int i = n-1; i>=0; i--){
            double suma1 = 0;
            cota[i] = cota[i] + Q[i].first;
            double den = 0;
            for(int r = 0; r<R; r++){
                paux[i][r] += sobra[r];
                sobra[r] = 0;
                double ganancia = W*log(s0[t][Q[i].second][bitmask[r]][r]*paux[i][r]);
                double perdida = 0;
                for(int j = 0; j<i; j++){
                    perdida -= W*d[Q[i].second][Q[j].second][bitmask[r]][r];
                    cota[j] -= W*d[Q[i].second][Q[j].second][bitmask[r]][r];
                }
                if(perdida < ganancia){
                    cota[i] = cota[i] + perdida;
                    suma1 = suma1 + ganancia;
                    den = den + 1;
                }
                else{
                    sobra[r] = paux[i][r];
                    paux[i][r] = 0;
                    for(int j = 0; j<i; j++){
                        cota[j] += W*d[Q[i].second][Q[j].second][bitmask[r]][r];
                    }
                }
            }

            if(suma1 >= cota[i]){
                sum++;
                double alpha = exp((cota[i] - suma1)/(W*den));
                for(int r = 0; r<R; r++){
                    if(paux[i][r] == 0) continue;
                    sobra[r] = (1.0 - alpha)*paux[i][r];
                    paux[i][r] = alpha*paux[i][r];
                }
            }
            else{
                fuera[i] = true;
                for(int r = 0; r<R; r++){
                    if(paux[i][r] == 0) continue;
                    sobra[r] = paux[i][r];
                    paux[i][r] = 0;
                    for(int j = 0; j<i; j++){
                        cota[j] += W*d[Q[i].second][Q[j].second][bitmask[r]][r];
                    }
                }
            }
        }

        if(sum > maxi){
            maxi = sum;
            bitmaskmax = bitmask;
            pauxmaxi = paux;
        }
    }

    if(maxi == 0) return {0, vector<vector<vector<double> > > (N, vector<vector<double> > (K, vector<double>(R, 0)))};

    vector<vector<vector<double> > > pfinal(N, vector<vector<double> > (K, vector<double>(R, 0)));
    
    for(int i = 0; i<n; i++){
        for(int r = 0; r<R; r++){
            pfinal[Q[i].second][bitmaskmax[r]][r] = pauxmaxi[i][r];
        }
    }

    return {maxi, pfinal};
}

pair<int, vector<vector<vector<double> > > >  Tunosolve2(vector<pair<double, int> > &Q, int t){
    
    double contar = 0;
    int n = Q.size();

    vector<bool> usado(R, false);
    vector<double> restante(K, R);
    vector<vector<double> > restante2(R, vector<double>(K, 4));
    vector<bool> listylor(n, false);
    vector<vector<vector<double> > > pfinal(N, vector<vector<double> > (K, vector<double>(R, 0)));
    vector<double> need(n);
    int libre = R;

    for(int prendidos = min(n, 3); prendidos >= 1; prendidos--){

        int bitmasknumero = 0;
        vector<int> bitmask(n, 0);
        for(int i = 0; i<prendidos; i++){
            bitmask[n-1-i] = 1;
        }
        do{
            vector<int> bits;
            bool saltar = false;
            for(int i = 0; i<n; i++){
                need[i] = 0;
                if(bitmask[i] == 1){
                    bits.push_back(i);
                    if(listylor[i] == true){
                        saltar = true;
                    }
                    bitmasknumero = bitmasknumero + (1<<i);
                }
            }
            if(saltar == true) continue;
            for(int r = 0; r<R; r++){
                if(usado[r] == true) continue;
                for(int k = 0; k<K; k++){
                    double sump = 0;
                    for(auto i: bits){
                        double sumae = 0;
                        for(auto j: bits){
                            if(j == i) continue;
                            sumae = sumae + d[Q[i].second][Q[j].second][k][r];
                        }
                        need[i] = expm1l(Q[i].first/(exp(sumae)*W))/s0[t][Q[i].second][k][r];
                        if(need[i] > 4){
                            sump = 30;
                        }
                        sump = sump + need[i];
                    }
                    if(sump < restante[k] && sump < restante2[r][k]){
                        libre--;
                        contar = contar + prendidos;
                        for(auto i: bits){
                            pfinal[Q[i].second][k][r] = need[i];
                            listylor[i] = true;
                            restante[k] -= pfinal[Q[i].second][k][r];
                            restante2[r][k] -= pfinal[Q[i].second][k][r];
                        }
                        usado[r] = true;
                        r = R;
                        k = K;
                        break;
                    }
                }
            } 
        }while(next_permutation(bitmask.begin(), bitmask.end()));
    }

    if(libre == 0) return {contar, pfinal};

    for(int r = 0; r<R; r++){
        if(usado[r] == true) continue;
        for(int i = 0; i<n; i++){
            if(listylor[i] == true) continue;
            vector<double> paux(K, 0);
            int contarK = 0;
            for(int k = 0; k<K; k++){
                paux[k] = min(restante2[r][k], restante[k]/double(libre));
                if(paux[k] != 0){
                    contarK++;
                }
            }
            if(contarK == 0) continue;
            double mean = 1;
            for(int k = 0; k<K; k++){
                if(paux[k] == 0) continue;
                mean = mean*pow(paux[k]*s0[t][Q[i].second][k][r], duno/double(contarK));
            }
            if(W*double(contarK)*log1p(mean) >= Q[i].first){
                contar++;
                listylor[i] = true;
                for(int k = 0; k<K; k++){
                    pfinal[Q[i].second][k][r] = paux[k];
                }
                usado[r] = true;
                break;
            }
        }
    }
    
    vector<vector<double> > gananciax(n, vector<double> (R, 0));
    vector<set<pair<double, int> > > ganancia(n);
    vector<double> librex(n, 0);
    vector<double> suma(n, 0);
    vector<double> cota(n);
    vector<int> quien(R, -1);
    
    for(int i = 0; i<n; i++){
        cota[i] = Q[i].first;
        if(listylor[i] == true) continue;
        for(int r = 0; r<R; r++){
            if(usado[r] == true) continue;
            vector<double> paux(K, 0);
            int contarK = 0;
            for(int k = 0; k<K; k++){
                paux[k] = min(restante2[r][k], restante[k]/double(libre));
                if(paux[k] != 0){
                    contarK++;
                }
            }
            if(contarK == 0) continue;
            double mean = 1;
            for(int k = 0; k<K; k++){
                if(paux[k] == 0) continue;
                mean = mean*pow(paux[k]*s0[t][Q[i].second][k][r], duno/double(contarK));
            }
            gananciax[i][r] = W*double(contarK)*log1p(mean);
            ganancia[i].insert({(-W*double(contarK)*log1p(mean))/Q[i].first, r});
            librex[i] = librex[i] + W*double(contarK)*log1p(mean);
        }
    }
    
    for(int i = 0; i<n; i++){

        while(ganancia[i].size() > 0){
            if(suma[i] >= cota[i]) break;
            
            auto nodo = *ganancia[i].begin();
            ganancia[i].erase(nodo);
            double profit = -nodo.first*Q[i].first;
            int r = nodo.second;
            int exid = quien[r];
            
            if(exid == -1){
                suma[i] = suma[i] + profit;
                quien[r] = i;
                for(int k = 0; k<K; k++){
                    pfinal[Q[i].second][k][r] = min(restante2[r][k], restante[k]/double(libre));
                }
                for(int ix = 0; ix<n; ix++){
                    if(listylor[ix] == true) continue;
                    librex[ix] = librex[ix] - gananciax[ix][r];
                }
                continue;
            }
        }
    
        if(suma[i] < cota[i]){
            for(int r = 0; r<R; r++){
                if(quien[r] == i){
                    quien[r] = -1;
                    for(int k = 0; k<K; k++){
                        pfinal[Q[i].second][k][r] = 0;
                    }
                }
            }
        }
        else{
            contar = contar + 1;
        }
    }
    
    return {contar, pfinal};
}

int Tuno(){

    for(int t = 0; t<T; t++){
        for(int k = 0; k<K; k++){
            for(int r = 0; r<R; r++){
                for(int i = 0; i<N; i++){
                    p[t][i][k][r] = 0;
                }
            }
        }
    }
    
    vector<vector<pair<double, int> > > Q(T);
    int contar = 0;
    for(auto u: A){
        Q[u[4] + rand()%u[5]].push_back({u[2], u[3]});
    }
    for(int t = 0; t<T; t++){
        if(Q[t].size() == 0) continue;
        pair<int, vector<vector<vector<double> > > > u = {0, vector<vector<vector<double> > > (0)};
        pair<int, vector<vector<vector<double> > > > v = {0, vector<vector<vector<double> > > (0)};
        u = Tunosolve(Q[t], t);
        if(Q[t].size() <= 18){
            v = Tunosolve2(Q[t], t);
        }
        if(u.first == 0 && v.first == 0){
            continue;
        }
        if(u.first >= v.first){
            contar = contar + u.first;
            for(int i = 0; i<N; i++){
                for(int k = 0; k<K; k++){
                    for(int r = 0; r<R; r++){
                        p[t][i][k][r] = u.second[i][k][r];
                    }
                }
            }
        }
        else{
            contar = contar + v.first;
            for(int i = 0; i<N; i++){
                for(int k = 0; k<K; k++){
                    for(int r = 0; r<R; r++){
                        p[t][i][k][r] = v.second[i][k][r];
                    }
                }
            }
        }
    }

    return contar;
}

int main(){

    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    string auxs;
    int maxT = 0;

    cin>>N>>K>>T>>R;

    s0.resize(T, vector<vector<vector<double> > >(N, vector<vector<double> > (K, vector<double>(R))));
    d.resize(N, vector<vector<vector<double> > >(N, vector<vector<double> > (K, vector<double>(R))));

    for(int t = 0; t<T; t++){
        for(int k = 0; k<K; k++){
            for(int r = 0; r<R; r++){
                for(int i = 0; i<N; i++){
                    cin>>auxs;
                    s0[t][i][k][r] = parseador(auxs);
                }
            }
        }
    }

    for(int k = 0; k<K; k++){
        for(int r = 0; r<R; r++){
            for(int j = 0; j<N; j++){
                for(int i = 0; i<N; i++){
                    cin>>auxs;
                    d[i][j][k][r] = parseador(auxs);
                }
            }
        }
    }

    cin>>J;
    A.resize(J, vector<int>(6));

    for(int i = 0; i<J; i++){
        cin>>A[i][1]>>A[i][2]>>A[i][3]>>A[i][4]>>A[i][5];
        maxT = max(maxT, A[i][5]);
    }

    p.resize(T, vector<vector<vector<double> > >(N, vector<vector<double> > (K, vector<double>(R))));

    vector<vector<vector<vector<double> > > > resp(T, vector<vector<vector<double> > >(N, vector<vector<double> > (K, vector<double>(R, 0))));
    double maxi = 0;
    double res = 0;
    
    res = LocalMaxSinInterUpdate(2);

    if(res > maxi){
        maxi = res;
        swap(resp, p);
    }

    if(maxT > 1){
        res = LocalMaxSinInterUpdate(5);
    }

    if(res > maxi){
        maxi = res;
        swap(resp, p);
    }

    if(maxT <= 3 && R <= 7){
        res = Tuno();
        if(res > maxi){
            maxi = res;
            swap(resp, p);
        }
    }

    cout<<fixed<<setprecision(14);

    for(int t = 0; t<T; t++){
        for(int k = 0; k<K; k++){
            for(int r = 0; r<R; r++){
                for(int i = 0; i<N; i++){
                    cout<<resp[t][i][k][r]<<" ";
                }
                cout<<"\n";
            }
        }
    }

    return 0;
}
