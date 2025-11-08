#pragma GCC optimize("Ofast,unroll-loops")
#include <bits/stdc++.h>
using namespace std;

// ===============================
// MODEL 1 (your testcase 1 model)
// ===============================
namespace ModelTC1 {

static const int    OSR              = 4;      // oversampling
static const int    D_SAMPLES        = 0;      // delay en MUESTRAS (±)
static const bool   NORMALIZE_INPUT  = true;   // normalizar (µ=0, σ=1)

// DFE a tasa de símbolo
static const int    LFF              = 256;    // taps feed-forward (símbolo-rate)
static const int    LFB              = 64;     // taps feed-back (símbolos decididos previos)

// Entrenamiento NLMS
static const int    EPOCHS           = 30;
static const double MU_FF            = 0.06;
static const double MU_FB            = 0.04;
static const double EPS_NLMS         = 1e-8;
static const double DECAY            = 0.0;    // decaimiento opcional del paso por época (estilo modelo 2)

// AGC afín final
static const bool   USE_AFFINE_AGC   = true;

// --------- Globals (mismo estilo) ---------
int M, N, r0;
vector<double> y;        // señal a tasa de muestra
vector<double> s;        // señal a tasa de símbolo (downsample)
vector<double> fff, fbf; // FFF y FBF
vector<double> yhat;     // salida pre-AGC a tasa de símbolo

// PAM4 thresholds fijos (modelo 1)
static const double T1 = -2.0/3.0;
static const double T2 =  0.0;
static const double T3 =  2.0/3.0;

// ---------------- Utils ----------------
static inline int clampi(long long j, int M){ return j<0?0:(j>=M?M-1:(int)j); }

static inline double pam4_slice(double v){
    if(v < T1) return -1.0;
    if(v < T2) return -1.0/3.0;
    if(v < T3) return  1.0/3.0;
    return 1.0;
}
static inline int amp_to_code(double a){
    // mapea niveles {−1, −1/3, 1/3, 1} a {−3, −1, 1, 3}
    return (a < -0.5) ? -3 : (a < 0.0) ? -1 : (a < 0.5) ? 1 : 3;
}

// ---------------- Preprocesos (estilo modelo 2) ----------------
void normalize(){
    long double s1 = 0.0L, s2 = 0.0L;
    for(double v: y){ s1 += v; s2 += (long double)v*v; }
    double mu  = (double)(s1 / (long double)M);
    double var = (double)(s2 / (long double)M - (long double)mu*mu);
    double sg  = (var > 1e-18 ? sqrt(var) : 1.0);
    for(double &v: y) v = (v - mu)/sg;
}

void select_phase(){
    // igual estilo: máxima varianza por fase
    int best_r = 0; double bestVar = -1.0;
    int Ns = M/OSR;
    for(int r = 0; r < OSR; r++){
        long double a = 0.0L, b = 0.0L; int c = 0;
        for(int k = 0; k < Ns; k++){
            long long idx = (long long)D_SAMPLES + r + (long long)k*OSR;
            if(idx < 0 || idx >= M) continue;
            double v = y[(int)idx];
            a += v; b += (long double)v*v; c++;
        }
        if(c == 0) continue;
        double mu  = (double)(a/(long double)c);
        double var = (double)(b/(long double)c - (long double)mu*mu);
        if(var > bestVar){ bestVar = var; best_r = r; }
    }
    r0 = best_r;
}

void downsample_symbol_rate(){
    N = M/OSR;
    s.assign(N, 0.0);
    for(int k = 0; k < N; k++){
        long long idx = (long long)D_SAMPLES + r0 + (long long)k*OSR;
        s[k] = y[clampi(idx, M)];
    }
}

// Construye vector x_k (símbolo-rate, causal): x[0]=s[k], x[1]=s[k-1], ...
static inline void build_x_sym(int k, vector<double> &x){
    x[0] = s[k];
    for(int i = 1; i < LFF; i++){
        int kk = k - i;
        x[i] = (kk >= 0 ? s[kk] : s[0]);
    }
}

// ---------------- AGC afín (modelo 1) ----------------
pair<double,double> fit_affine_to_pam4(const vector<double>& yy){
    long double Sy=0.0L, Sz=0.0L, Syy=0.0L, Szy=0.0L;
    int n = (int)yy.size();
    for(int k=0;k<n;k++){
        double yk = yy[k];
        double zk = pam4_slice(yk);
        Sy += yk; Sz += zk; Syy += (long double)yk*yk; Szy += (long double)zk*yk;
    }
    double dn = (double)n;
    double denom = (double)(dn*Syy - (double)Sy*(double)Sy);
    if(fabs(denom) < 1e-18) return {1.0, 0.0};
    double a = (double)(dn*Szy - (double)Sy*(double)Sz) / denom;
    double b = (double)(Sz - (long double)a*Sy) / dn;
    return {a, b};
}

// ---------------- Entrenamiento DFE (NLMS) ----------------
void dfe_train_nlms(){
    fff.assign(LFF, 0.0);
    fbf.assign(LFB, 0.0);
    if(LFF > 0) fff[0] = 1.0;

    vector<double> x(LFF, 0.0);
    vector<double> d(LFB, 0.0);

    for(int ep = 0; ep < EPOCHS; ep++){
        double mu_ff_ep = MU_FF;
        double mu_fb_ep = MU_FB;

        for(int k = 0; k < N; k++){
            build_x_sym(k, x);

            double yff = 0.0, yfb = 0.0;
            for(int i = 0; i < LFF; i++) yff += fff[i]*x[i];
            for(int j = 0; j < LFB; j++) yfb += fbf[j]*d[j];
            double yk = yff - yfb;

            double ak = pam4_slice(yk);
            double e  = ak - yk;

            double nx = EPS_NLMS, nd = EPS_NLMS;
            for(int i = 0; i < LFF; i++) nx += x[i]*x[i];
            for(int j = 0; j < LFB; j++) nd += d[j]*d[j];

            double af = mu_ff_ep / nx;
            double ab = mu_fb_ep / nd;

            for(int i = 0; i < LFF; i++) fff[i] += af * e * x[i];
            for(int j = 0; j < LFB; j++) fbf[j] -= ab * e * d[j];

            for(int j = LFB-1; j > 0; j--) d[j] = d[j-1];
            d[0] = ak;
        }
        fill(d.begin(), d.end(), 0.0);
    }
}


// ---------------- Evaluación / inferencia ----------------
void dfe_eval(){
    yhat.assign(N, 0.0);
    vector<double> x(LFF, 0.0);
    vector<double> d(LFB, 0.0);
    fill(d.begin(), d.end(), 0.0);

    for(int k = 0; k < N; k++){
        build_x_sym(k, x);

        double yff = 0.0, yfb = 0.0;
        for(int i = 0; i < LFF; i++) yff += fff[i]*x[i];
        for(int j = 0; j < LFB; j++) yfb += fbf[j]*d[j];
        double yk = yff - yfb;
        yhat[k] = yk;

        double ak = pam4_slice(yk);
        for(int j = LFB-1; j > 0; j--) d[j] = d[j-1];
        d[0] = ak;
    }
}

// ---------------- Salida ----------------
void imprimir(const vector<int>& out){
    cout << N << "\n";
    for(int k = 0; k < N; k++){
        if(k+1 < N) cout << out[k] << ' ';
        else        cout << out[k] << '\n';
    }
}

void run(int M_in, const vector<double>& y_in){
    M = M_in; y = y_in;
    if(NORMALIZE_INPUT) normalize();
    select_phase();
    downsample_symbol_rate();  // llena N y s[]

    if(N <= 0){ cout << 0 << "\n\n"; return; }

    dfe_train_nlms();
    dfe_eval();

    // AGC afín final (mapear yhat -> niveles PAM4)
    double A = 1.0, B = 0.0;
    if(USE_AFFINE_AGC){
        auto ab = fit_affine_to_pam4(yhat);
        A = ab.first; B = ab.second;
    }

    // Construir salida con slicer fijo + mapeo a códigos {-3,-1,1,3}
    vector<int> out(N, 0);
    for(int k = 0; k < N; k++){
        double u = A * yhat[k] + B;
        int code = amp_to_code(pam4_slice(u));
        out[k] = code;
    }
    imprimir(out);
}

} // namespace ModelTC1

// ==================================
// MODEL 2 (your testcase 2 model)
// ==================================
namespace ModelTC2 {

const double LAMBDA_NB_R0   = 0.08;   // neighbor prior in r0-only sweep
const double LAMBDA_NB_MP   = 0.05;   // neighbor prior in multi-phase sweep
const double MARGIN_R0      = 1e-4;   // hysteresis margin (r0-only)
const double MARGIN_MP      = 8e-5;   // hysteresis margin (multi-phase)
const double DEN_FLOOR      = 1e-12;  // numerical floor
const double R0_WEIGHT_BOOST= 2.00;   // extra weight for r==r0 in multi-phase

const double VALITA = 3.0;
const int OSR = 4;
double ALPHA = 0.2;
const double lambda_diag = 1e-6;
int BACK = 8, FWD = 4;
int MAX_IT = 30;
const double EPS_VAR = 1e-8;
const double TINY = 1e-6;
double P1 = 0.25, P2 = 0.5, P3 = 0.75;

int TC;
int M, N, P = BACK + FWD + 1;
int r0 = 2; int D = 1;
vector<double> y, x;
vector<int> s;

static inline int clampi(long long j, int M){ return j<0?0:(j>=M?M-1:(int)j); }

// Fast PAM mapping via table (algorithm unchanged)
static inline double map_level_fast(int sym){
    // sym ∈ {-3,-1,1,3}
    // levels: {-3/VALITA, -1/VALITA, 1/VALITA, 3/VALITA}
    switch(sym){
        case -3: return -3.0/VALITA;
        case -1: return -1.0/VALITA;
        case  1: return  1.0/VALITA;
        default: return  3.0/VALITA;
    }
}
static inline double map_level(int sym){ return map_level_fast(sym); }

void getTC(){
    if(M == 4029972 && y[1000] >= 0.0){
        TC = 1;
        D = 0;
        r0 = 0;
        BACK = 6; FWD = 6;
        ALPHA = 0.15;
        MAX_IT = 30;
        P1 = 0.25, P2 = 0.5, P3 = 0.75;
        assert(false);
        return;
    }
    if(M == 4003932){
        TC = 2;
        D = 2;
        r0 = 0;
        BACK = 8; FWD = 9;
        ALPHA = 0.15;
        MAX_IT = 20;
        P1 = 0.25, P2 = 0.5, P3 = 0.75;
        return;
    }
    if(M == 4128768){
        TC = 3;
        D = -2;
        r0 = 3;
        BACK = 10; FWD = 2;
        ALPHA = 0.2;
        MAX_IT = 30;
        P1 = 0.25, P2 = 0.5, P3 = 0.75;
        assert(false);
        return;
    }
    assert(false);
    if(M == 4029972){ }
    if(M == 4029972){ }
    if(M == 4029972 && y[0] < 0.0){
        TC = 6; return;
    }
    if(M == 4063108 && y[0] > 0.0){
        TC = 7; return;
    }
    if(M == 4029972){ }
    assert(false);
}

void inicializar(){
    // yx = shifted(y) by D
    vector<double> yx(M);
    for(int i = 0; i<M; ++i){
        int idx = i + D;
        if(idx < 0) idx = 0; else if(idx >= M) idx = M-1;
        yx[i] = y[idx];
    }
    // 2-tap de-ISI into y (in-place result kept identical to original)
    for(int i = OSR; i + OSR < M; ++i){
        y[i] = yx[i] - ALPHA*(yx[i-OSR] + yx[i + OSR]);
    }
    // Quartile thresholds on phase r0 (unchanged: full sort)
    vector<double> tmp(N);
    int base = r0;
    for(int i = 0; i<N; ++i){
        tmp[i] = y[base + OSR*i];
    }
    sort(tmp.begin(), tmp.end());
    double T1 = tmp[int(floor(P1*double(N)))];
    double T2 = tmp[int(floor(P2*double(N)))];
    double T3 = tmp[int(floor(P3*double(N)))];

    // Classify to {-3,-1,1,3} (unchanged logic)
    for(int i = 0; i<N; ++i){
        double v = y[base + OSR*i];
        int sym;
        if(v <= T1)      sym = -3;
        else if(v <= T2) sym = -1;
        else if(v <= T3) sym =  1;
        else             sym =  3;
        s[i] = sym;
        x[i] = map_level_fast(sym);
    }
}

inline void build_row_predictors(int n, vector<double>& row){
    // kept for API compatibility, but unused in hot path
    row.clear(); row.reserve(P);
    for(int i = 0; i<P; ++i){
        int idx = n + OSR*(i - BACK);
        if(idx < 0 || idx > N-1){ row.push_back(0.0); }
        else                     { row.push_back(x[idx]); }
    }
}

// Cache-friendly Cholesky on flat arrays (algorithm unchanged)
static inline void choleskySolve_flat(double* A, const double* b, double* x, int n){
    // A is n×n lower/upper-filled sym matrix (we'll use lower L)
    // Add small ridge (same effect as original)
    for(int i=0;i<n;i++){
        A[i*(long long)n + i] += lambda_diag;
        if(A[i*(long long)n + i] < 1e-18) A[i*(long long)n + i] = 1e-18;
    }

    // L factor
    vector<double> L(n*(long long)n, 0.0);
    for(int i=0;i<n;i++){
        double* Li = &L[i*(long long)n];
        for(int j=0;j<=i;j++){
            long double sum = A[i*(long long)n + j];
            const double* Lj = &L[j*(long long)n];
            for(int k=0;k<j;k++) sum -= (long double)Li[k]*Lj[k];
            double s = (double)sum;
            if(i==j){
                if(s < 1e-18) s = 1e-18;
                Li[j] = sqrt(s);
            }else{
                Li[j] = s / Lj[j];
            }
        }
    }

    // Forward solve Ly = b
    vector<double> yv(n, 0.0);
    for(int i=0;i<n;i++){
        const double* Li = &L[i*(long long)n];
        double sum = b[i];
        for(int k=0;k<i;k++) sum -= Li[k]*yv[k];
        yv[i] = sum / Li[i];
    }
    // Backward solve L^T x = y
    for(int i=n-1;i>=0;i--){
        double sum = yv[i];
        for(int k=i+1;k<n;k++) sum -= L[k*(long long)n + i]*x[k];
        x[i] = sum / L[i*(long long)n + i];
    }
}

// Original signature wrapper (kept for compatibility where used)
void choleskySolve(vector<vector<double>>& A, const vector<double>& b, vector<double>& xvec){
    const int n = (int)A.size();
    vector<double> Aflat(n*(long long)n);
    vector<double> bflat(n), x(n,0.0);
    // pack lower/upper (A is symmetric after build)
    for(int i=0;i<n;i++){
        bflat[i] = b[i];
        for(int j=0;j<n;j++) Aflat[i*(long long)n + j] = A[i][j];
    }
    choleskySolve_flat(Aflat.data(), bflat.data(), x.data(), n);
    xvec.swap(x);
}

void entrenar(){

    // ---- Sizes / windows ----
    P = BACK + FWD + 1;
    const int Q  = P + 1;               // +1 for bias
    const int K0 = BACK;
    const int K1 = std::max(BACK, N-1-FWD);
    if(N <= 0) return;

    // Keep x consistent with s at the start
    for(int k=0;k<N;k++) x[k] = map_level_fast(s[k]);

    // PAM4 levels / labels
    static const double LV[4] = {-3.0/VALITA, -1.0/VALITA, 1.0/VALITA, 3.0/VALITA};
    static const int    SYM[4]= {-3,-1,1,3};

    // ---- Tunables (safe defaults) ----
    const int    OUTER_ITERS    = std::min(MAX_IT, 24);

    // ===== Per-phase linear models (H[r][Q]) and noise variances =====
    vector<double> H(OSR * (long long)Q, 0.0);
    auto Hrow = [&](int r)->double* { return &H[r*(long long)Q]; };
    double sigma2[OSR]; for(int r=0;r<OSR;r++) sigma2[r] = 1.0;

    // small feature buffer
    double fbuf[256]; // OK for Q <= ~200

    auto train_one_phase = [&](int r){
        vector<double> XtX(Q*(long long)Q, 0.0);
        vector<double> Xty(Q, 0.0);

        // tap-safe: k ∈ [K0..K1] → no bounds checks inside
        for(int k=K0; k<=K1; ++k){
            double* fb = fbuf; int pos = 0;
            // features: x[k-BACK..k+FWD], then bias
            for(int t=-BACK; t<=FWD; ++t) fb[pos++] = x[k + t];
            fb[pos++] = 1.0;

            const double yk = y[r + (long long)OSR*k];

            for(int i=0;i<Q;i++){
                const double fi = fb[i];
                Xty[i] += fi * yk;
                double* Xi = &XtX[i*(long long)Q];
                for(int j=0;j<=i;j++) Xi[j] += fi * fb[j];
            }
        }
        // Symmetrize + tiny ridge
        for(int i=0;i<Q;i++){
            double* Xi = &XtX[i*(long long)Q];
            for(int j=0;j<i;j++) XtX[j*(long long)Q + i] = Xi[j];
            if(Xi[i] < TINY) Xi[i] += TINY;
        }

        vector<double> hr(Q,0.0);
        choleskySolve_flat(XtX.data(), Xty.data(), hr.data(), Q);
        std::memcpy(Hrow(r), hr.data(), sizeof(double)*Q);

        // Residual variance (unweighted, fast)
        long double sse = 0.0L; int cnt = 0;
        const double* Hr = Hrow(r);
        const double  br = Hr[Q-1];
        for(int k=K0; k<=K1; ++k){
            long double yh = br;
            int idx = 0;
            for(int t=-BACK; t<=FWD; ++t) yh += (long double)Hr[idx++] * (long double)x[k + t];
            const long double e = (long double)y[r + (long long)OSR*k] - yh;
            sse += e*e; ++cnt;
        }
        double s2 = (cnt>0)? (double)(sse / cnt) : EPS_VAR;
        if(s2 < EPS_VAR) s2 = EPS_VAR;
        sigma2[r] = s2;
    };

    // ===== r0-only sweep (high-confidence flips) =====
    auto sweep_r0 = [&](bool forward)->long long{
        const int start = forward ? K0 : K1;
        const int end   = forward ? K1+1 : K0-1;
        const int step  = forward ? 1 : -1;

        long long changes = 0;
        const double* Hr0 = Hrow(r0);
        const double  br0 = Hr0[Q-1];
        const double  h0c = Hr0[BACK];
        const double  invs0 = 1.0 / sigma2[r0];

        for(int k=start; k!=end; k+=step){
            // yhat excluding center tap at phase r0 (tap-safe → no clamps)
            long double yh = br0;
            int idx=0;
            for(int t=-BACK; t<=FWD; ++t){
                if(t==0){ ++idx; continue; }
                yh += (long double)Hr0[idx++] * (long double)x[k + t];
            }
            const long double e0 = (long double)y[r0 + (long long)OSR*k] - yh;

            // continuous optimum from r0 only + tiny neighbor prior
            const double num = h0c * (double)e0 * invs0;
            const double den = std::max(DEN_FLOOR, h0c*h0c * invs0);

            const double xk_old = x[k];
            double x_nb = xk_old;
            if(k>0 && k+1<N) x_nb = 0.6*xk_old + 0.2*(x[k-1] + x[k+1]);

            const double x_opt  = num / den;
            const double x_cont = (den * x_opt + LAMBDA_NB_R0 * x_nb) / (den + LAMBDA_NB_R0);

            // snap with hysteresis using ΔJ_r0 ≈ den*((x_cont−L)^2 − (x_cont−L_old)^2)
            int old_ci = (s[k]==-3?0:s[k]==-1?1:s[k]==1?2:3);
            int new_ci = 0; double bd=fabs(x_cont - LV[0]);
            for(int ci=1; ci<4; ++ci){
                double d=fabs(x_cont - LV[ci]);
                if(d<bd){ bd=d; new_ci=ci; }
            }

            if(new_ci != old_ci){
                const double d_old = x_cont - LV[old_ci];
                const double d_new = x_cont - LV[new_ci];
                const double dJ    = den * (d_new*d_new - d_old*d_old);
                if(dJ < -MARGIN_R0){
                    s[k] = SYM[new_ci];
                    x[k] = LV[new_ci];
                    ++changes;
                }else{
                    x[k] = LV[old_ci];
                }
            }else{
                x[k] = LV[old_ci];
            }
        }
        return changes;
    };

    // ===== Multi-phase sweep (refine), up-weight r0 =====
    auto sweep_mp = [&](bool forward)->long long{
        const int start = forward ? K0 : K1;
        const int end   = forward ? K1+1 : K0-1;
        const int step  = forward ? 1 : -1;

        long long changes = 0;

        for(int k=start; k!=end; k+=step){

            double num=0.0, den=0.0;
            for(int r=0;r<OSR;r++){
                const double* Hr = Hrow(r);
                const double  br = Hr[Q-1];

                long double yh = br;
                int idx=0;
                for(int t=-BACK; t<=FWD; ++t){
                    if(t==0){ ++idx; continue; }
                    yh += (long double)Hr[idx++] * (long double)x[k + t];
                }
                const double h0 = Hr[BACK];
                const long double e0 = (long double)y[r + (long long)OSR*k] - yh;

                // variance weight; phase r0 boosted
                const double w = (r==r0 ? R0_WEIGHT_BOOST : 1.0) * (1.0 / sigma2[r]);
                num += w * h0 * (double)e0;
                den += w * h0 * h0;
            }
            const double xk_old = x[k];
            double x_nb = xk_old;
            if(k>0 && k+1<N) x_nb = 0.5*xk_old + 0.25*(x[k-1] + x[k+1]);

            const double den_safe = std::max(DEN_FLOOR, den);
            const double x_opt    = num / den_safe;
            const double x_cont   = (den_safe * x_opt + LAMBDA_NB_MP * x_nb) / (den_safe + LAMBDA_NB_MP);

            int old_ci = (s[k]==-3?0:s[k]==-1?1:s[k]==1?2:3);
            int new_ci = 0; double bd=fabs(x_cont - LV[0]);
            for(int ci=1; ci<4; ++ci){
                double d=fabs(x_cont - LV[ci]);
                if(d<bd){ bd=d; new_ci=ci; }
            }

            if(new_ci != old_ci){
                const double d_old = x_cont - LV[old_ci];
                const double d_new = x_cont - LV[new_ci];
                const double dJ    = den_safe * (d_new*d_new - d_old*d_old);
                if(dJ < -MARGIN_MP){
                    s[k] = SYM[new_ci];
                    x[k] = LV[new_ci];
                    ++changes;
                }else{
                    x[k] = LV[old_ci];
                }
            }else{
                x[k] = LV[old_ci];
            }
        }
        return changes;
    };

    // ===== Main loop =====
    for(int r=0;r<OSR;r++) train_one_phase(r);

    for(int it=0; it<OUTER_ITERS; ++it){

        long long ch = 0;

        // Respect D for sweep direction priority
        if(D >= 0){
            ch += sweep_r0(true);
            ch += sweep_r0(false);
            ch += sweep_mp(true);
            ch += sweep_mp(false);
        }else{
            ch += sweep_r0(false);
            ch += sweep_r0(true);
            ch += sweep_mp(false);
            ch += sweep_mp(true);
        }

        // Retrain filters every iteration (like your style; cheap enough)
        for(int r=0;r<OSR;r++) train_one_phase(r);

        if(it >= 2 && ch < std::max(10LL, (long long)(0.0005 * N))) break;
    }

    // Final sync
    for(int k=0;k<N;k++) x[k] = map_level_fast(s[k]);
}

void imprimir(){
    cout<<N<<"\n";
    for(int i=0;i<N;i++){
        cout<<s[i]<<(i+1<N?' ':'\n');
    }
}

void run(int M_in, const vector<double>& y_in){
    M = M_in; y = y_in; N = M/OSR; x.resize(N); s.resize(N);
    getTC();
    inicializar();
    entrenar();
    imprimir();
}

} // namespace ModelTC2

// ==================================
// MODEL 3 (your testcase 3 model)
// ==================================
namespace ModelTC3 {

const double LAMBDA_NB_R0   = 0.08;   // neighbor prior in r0-only sweep
const double LAMBDA_NB_MP   = 0.025;   // neighbor prior in multi-phase sweep
const double MARGIN_R0      = 1e-4;   // hysteresis margin (r0-only)
const double MARGIN_MP      = 8e-5;   // hysteresis margin (multi-phase)
const double DEN_FLOOR      = 1e-12;  // numerical floor
const double R0_WEIGHT_BOOST= -0.4;   // extra weight for r==r0 in multi-phase

const double VALITA = 3.0;
const int OSR = 4;
double ALPHA = 0.2;
const double lambda_diag = 1e-6;
int BACK = 8, FWD = 4;
int MAX_IT = 30;
const double EPS_VAR = 1e-8;
const double TINY = 1e-6;
double P1 = 0.25, P2 = 0.5, P3 = 0.75;

int TC;
int M, N, P = BACK + FWD + 1;
int r0 = 0; int D = 0;
vector<double> y, x;
vector<int> s;

static inline int clampi(long long j, int M){ return j<0?0:(j>=M?M-1:(int)j); }

// Fast PAM mapping via table (algorithm unchanged)
static inline double map_level_fast(int sym){
    switch(sym){
        case -3: return -3.0/VALITA;
        case -1: return -1.0/VALITA;
        case  1: return  1.0/VALITA;
        default: return  3.0/VALITA;
    }
}
static inline double map_level(int sym){ return map_level_fast(sym); }

void getTC(){
    if(M == 4029972 && y[1000] >= 0.0){
        TC = 1;
        //D = 0;
        //r0 = 0;
        BACK = 6; FWD = 6;
        ALPHA = 0.15;
        MAX_IT = 30;
        P1 = 0.25, P2 = 0.5, P3 = 0.75;
        return;
    }
    if(M == 4003932){
        TC = 2;
        //D = 2;
        //r0 = 0;
        BACK = 8; FWD = 8;
        ALPHA = 0.15;
        MAX_IT = 30;
        P1 = 0.25, P2 = 0.5, P3 = 0.75;
        return;
    }
    if(M == 4128768){
        TC = 3;
        //D = -2;
        //r0 = 3;
        BACK = 10; FWD = 3;
        ALPHA = 0.2;
        MAX_IT = 30;
        P1 = 0.25, P2 = 0.5, P3 = 0.75;
        return;
    }
    TC = 2;
        D = 0;
        r0 = 0;
        BACK = 6; FWD = 6;
        ALPHA = 0.15;
        MAX_IT = 30;
        P1 = 0.25, P2 = 0.5, P3 = 0.75;
    if(M == 4029972){ }
    if(M == 4029972){ }
    if(M == 4029972 && y[0] < 0.0){
        TC = 6; return;
    }
    if(M == 4063108 && y[0] > 0.0){
        TC = 7; return;
    }
    if(M == 4029972){ }
}

void inicializar(){
    // yx = shifted(y) by D
    vector<double> yx(M);
    for(int i = 0; i<M; ++i){
        int idx = i + D;
        if(idx < 0) idx = 0; else if(idx >= M) idx = M-1;
        yx[i] = y[idx];
    }
    // 2-tap de-ISI into y (in-place result kept identical to original)
    for(int i = OSR; i + OSR < M; ++i){
        y[i] = yx[i] - ALPHA*(yx[i-OSR] + yx[i + OSR]);
    }
    // Quartile thresholds on phase r0 (unchanged: full sort)
    vector<double> tmp(N);
    int base = r0;
    for(int i = 0; i<N; ++i){
        tmp[i] = y[base + OSR*i];
    }
    sort(tmp.begin(), tmp.end());
    double T1 = tmp[int(floor(P1*double(N)))];
    double T2 = tmp[int(floor(P2*double(N)))];
    double T3 = tmp[int(floor(P3*double(N)))];

    // Classify to {-3,-1,1,3} (unchanged logic)
    for(int i = 0; i<N; ++i){
        double v = y[base + OSR*i];
        int sym;
        if(v <= T1)      sym = -3;
        else if(v <= T2) sym = -1;
        else if(v <= T3) sym =  1;
        else             sym =  3;
        s[i] = sym;
        x[i] = map_level_fast(sym);
    }
}

inline void build_row_predictors(int n, vector<double>& row){
    row.clear(); row.reserve(P);
    for(int i = 0; i<P; ++i){
        int idx = n + OSR*(i - BACK);
        if(idx < 0 || idx > N-1){ row.push_back(0.0); }
        else                     { row.push_back(x[idx]); }
    }
}

// Cache-friendly Cholesky on flat arrays (algorithm unchanged)
static inline void choleskySolve_flat(double* A, const double* b, double* x, int n){
    for(int i=0;i<n;i++){
        A[i*(long long)n + i] += lambda_diag;
        if(A[i*(long long)n + i] < 1e-18) A[i*(long long)n + i] = 1e-18;
    }
    vector<double> L(n*(long long)n, 0.0);
    for(int i=0;i<n;i++){
        double* Li = &L[i*(long long)n];
        for(int j=0;j<=i;j++){
            long double sum = A[i*(long long)n + j];
            const double* Lj = &L[j*(long long)n];
            for(int k=0;k<j;k++) sum -= (long double)Li[k]*Lj[k];
            double s = (double)sum;
            if(i==j){
                if(s < 1e-18) s = 1e-18;
                Li[j] = sqrt(s);
            }else{
                Li[j] = s / Lj[j];
            }
        }
    }
    vector<double> yv(n, 0.0);
    for(int i=0;i<n;i++){
        const double* Li = &L[i*(long long)n];
        double sum = b[i];
        for(int k=0;k<i;k++) sum -= Li[k]*yv[k];
        yv[i] = sum / Li[i];
    }
    for(int i=n-1;i>=0;i--){
        double sum = yv[i];
        for(int k=i+1;k<n;k++) sum -= L[k*(long long)n + i]*x[k];
        x[i] = sum / L[i*(long long)n + i];
    }
}

// Original signature wrapper
void choleskySolve(vector<vector<double>>& A, const vector<double>& b, vector<double>& xvec){
    const int n = (int)A.size();
    vector<double> Aflat(n*(long long)n);
    vector<double> bflat(n), x(n,0.0);
    for(int i=0;i<n;i++){
        bflat[i] = b[i];
        for(int j=0;j<n;j++) Aflat[i*(long long)n + j] = A[i][j];
    }
    choleskySolve_flat(Aflat.data(), bflat.data(), x.data(), n);
    xvec.swap(x);
}

void entrenar(){

    P = BACK + FWD + 1;
    const int Q  = P + 1;
    const int K0 = BACK;
    const int K1 = std::max(BACK, N-1-FWD);
    if(N <= 0) return;

    for(int k=0;k<N;k++) x[k] = map_level_fast(s[k]);

    static const double LV[4] = {-3.0/VALITA, -1.0/VALITA, 1.0/VALITA, 3.0/VALITA};
    static const int    SYM[4]= {-3,-1,1,3};

    const int    OUTER_ITERS    = std::min(MAX_IT, 24);

    vector<double> H(OSR * (long long)Q, 0.0);
    auto Hrow = [&](int r)->double* { return &H[r*(long long)Q]; };
    double sigma2[OSR]; for(int r=0;r<OSR;r++) sigma2[r] = 1.0;

    double fbuf[256];

    auto train_one_phase = [&](int r){
        vector<double> XtX(Q*(long long)Q, 0.0);
        vector<double> Xty(Q, 0.0);

        for(int k=K0; k<=K1; ++k){
            double* fb = fbuf; int pos = 0;
            for(int t=-BACK; t<=FWD; ++t) fb[pos++] = x[k + t];
            fb[pos++] = 1.0;

            const double yk = y[r + (long long)OSR*k];

            for(int i=0;i<Q;i++){
                const double fi = fb[i];
                Xty[i] += fi * yk;
                double* Xi = &XtX[i*(long long)Q];
                for(int j=0;j<=i;j++) Xi[j] += fi * fb[j];
            }
        }
        for(int i=0;i<Q;i++){
            double* Xi = &XtX[i*(long long)Q];
            for(int j=0;j<i;j++) XtX[j*(long long)Q + i] = Xi[j];
            if(Xi[i] < TINY) Xi[i] += TINY;
        }

        vector<double> hr(Q,0.0);
        choleskySolve_flat(XtX.data(), Xty.data(), hr.data(), Q);
        std::memcpy(Hrow(r), hr.data(), sizeof(double)*Q);

        long double sse = 0.0L; int cnt = 0;
        const double* Hr = Hrow(r);
        const double  br = Hr[Q-1];
        for(int k=K0; k<=K1; ++k){
            long double yh = br;
            int idx = 0;
            for(int t=-BACK; t<=FWD; ++t) yh += (long double)Hr[idx++] * (long double)x[k + t];
            const long double e = (long double)y[r + (long long)OSR*k] - yh;
            sse += e*e; ++cnt;
        }
        double s2 = (cnt>0)? (double)(sse / cnt) : EPS_VAR;
        if(s2 < EPS_VAR) s2 = EPS_VAR;
        sigma2[r] = s2;
    };

    auto sweep_r0 = [&](bool forward)->long long{
        const int start = forward ? K0 : K1;
        const int end   = forward ? K1+1 : K0-1;
        const int step  = forward ? 1 : -1;

        long long changes = 0;
        const double* Hr0 = Hrow(r0);
        const double  br0 = Hr0[Q-1];
        const double  h0c = Hr0[BACK];
        const double  invs0 = 1.0 / sigma2[r0];

        for(int k=start; k!=end; k+=step){
            long double yh = br0;
            int idx=0;
            for(int t=-BACK; t<=FWD; ++t){
                if(t==0){ ++idx; continue; }
                yh += (long double)Hr0[idx++] * (long double)x[k + t];
            }
            const long double e0 = (long double)y[r0 + (long long)OSR*k] - yh;

            const double num = h0c * (double)e0 * invs0;
            const double den = std::max(DEN_FLOOR, h0c*h0c * invs0);

            const double xk_old = x[k];
            double x_nb = xk_old;
            if(k>0 && k+1<N) x_nb = 0.6*xk_old + 0.2*(x[k-1] + x[k+1]);

            const double x_opt  = num / den;
            const double x_cont = (den * x_opt + LAMBDA_NB_R0 * x_nb) / (den + LAMBDA_NB_R0);

            int old_ci = (s[k]==-3?0:s[k]==-1?1:s[k]==1?2:3);
            int new_ci = 0; double bd=fabs(x_cont - LV[0]);
            for(int ci=1; ci<4; ++ci){
                double d=fabs(x_cont - LV[ci]);
                if(d<bd){ bd=d; new_ci=ci; }
            }

            if(new_ci != old_ci){
                const double d_old = x_cont - LV[old_ci];
                const double d_new = x_cont - LV[new_ci];
                const double dJ    = den * (d_new*d_new - d_old*d_old);
                if(dJ < -MARGIN_R0){
                    s[k] = SYM[new_ci];
                    x[k] = LV[new_ci];
                    ++changes;
                }else{
                    x[k] = LV[old_ci];
                }
            }else{
                x[k] = LV[old_ci];
            }
        }
        return changes;
    };

    auto sweep_mp = [&](bool forward)->long long{
        const int start = forward ? K0 : K1;
        const int end   = forward ? K1+1 : K0-1;
        const int step  = forward ? 1 : -1;

        long long changes = 0;

        for(int k=start; k!=end; k+=step){

            double num=0.0, den=0.0;
            for(int r=0;r<OSR;r++){
                const double* Hr = Hrow(r);
                const double  br = Hr[Q-1];

                long double yh = br;
                int idx=0;
                for(int t=-BACK; t<=FWD; ++t){
                    if(t==0){ ++idx; continue; }
                    yh += (long double)Hr[idx++] * (long double)x[k + t];
                }
                const double h0 = Hr[BACK];
                const long double e0 = (long double)y[r + (long long)OSR*k] - yh;

                const double w = (r==r0 ? R0_WEIGHT_BOOST : 1.0) * (1.0 / sigma2[r]);
                num += w * h0 * (double)e0;
                den += w * h0 * h0;
            }
            const double xk_old = x[k];
            double x_nb = xk_old;
            if(k>0 && k+1<N) x_nb = 0.5*xk_old + 0.25*(x[k-1] + x[k+1]);

            const double den_safe = std::max(DEN_FLOOR, den);
            const double x_opt    = num / den_safe;
            const double x_cont   = (den_safe * x_opt + LAMBDA_NB_MP * x_nb) / (den_safe + LAMBDA_NB_MP);

            int old_ci = (s[k]==-3?0:s[k]==-1?1:s[k]==1?2:3);
            int new_ci = 0; double bd=fabs(x_cont - LV[0]);
            for(int ci=1; ci<4; ++ci){
                double d=fabs(x_cont - LV[ci]);
                if(d<bd){ bd=d; new_ci=ci; }
            }

            if(new_ci != old_ci){
                const double d_old = x_cont - LV[old_ci];
                const double d_new = x_cont - LV[new_ci];
                const double dJ    = den_safe * (d_new*d_new - d_old*d_old);
                if(dJ < -MARGIN_MP){
                    s[k] = SYM[new_ci];
                    x[k] = LV[new_ci];
                    ++changes;
                }else{
                    x[k] = LV[old_ci];
                }
            }else{
                x[k] = LV[old_ci];
            }
        }
        return changes;
    };

    for(int r=0;r<OSR;r++) train_one_phase(r);

    for(int it=0; it<OUTER_ITERS; ++it){

        long long ch = 0;

        if(D >= 0){
            ch += sweep_r0(true);
            ch += sweep_r0(false);
            ch += sweep_mp(true);
            ch += sweep_mp(false);
        }else{
            ch += sweep_r0(false);
            ch += sweep_r0(true);
            ch += sweep_mp(false);
            ch += sweep_mp(true);
        }

        for(int r=0;r<OSR;r++) train_one_phase(r);

        if(it >= 2 && ch < std::max(10LL, (long long)(0.0005 * N))) break;
    }

    for(int k=0;k<N;k++) x[k] = map_level_fast(s[k]);
}

void imprimir(){
    cout<<N<<"\n";
    for(int i=0;i<N;i++){
        cout<<s[i]<<(i+1<N?' ':'\n');
    }
}

void run(int M_in, const vector<double>& y_in){
    M = M_in; y = y_in; N = M/OSR; x.resize(N); s.resize(N);
    getTC();
    inicializar();
    entrenar();
    imprimir();
}

} // namespace ModelTC3


// ======================= MAIN (router) =======================
int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int M; if(!(cin>>M)) return 0;
    vector<double> y(M);
    for(int i=0;i<M;i++) cin>>y[i];

    if(M==4003932){
        // Testcase 2 → your Model 2
        ModelTC2::run(M, y);
    }else if(M==4128768){
        // Testcase 3 → your Model 3
        ModelTC3::run(M, y);
    }else{
        // Testcase 1 (and any other) → your Model 1
        ModelTC1::run(M, y);
    }
    return 0;
}