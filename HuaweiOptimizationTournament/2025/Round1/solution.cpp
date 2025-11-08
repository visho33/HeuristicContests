#pragma GCC optimize("Ofast")
#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
typedef pair<ll, ll> ii;

/*
Residual edge used by Dinic.
- v: tail node (from) of this directed residual edge
- u: head node (to) of this directed residual edge
- cap: capacity remaining on this edge
- flow: current flow sent on this edge (cap - flow is residual)
- cost: per-unit cost (used to prefer cheaper augmenting paths)
*/
struct FlowEdge{
    int v, u;
    long long cap, flow;
    long long cost;
    FlowEdge(int v,int u,long long cap,long long cost=0)
        : v(v), u(u), cap(cap), flow(0), cost(cost) {}
};

/*
Max-flow engine with epsilon-scheduled reduced-cost filtering to bias augmenting paths toward cheaper arcs.
- edges: list of residual edges (paired forward/backward)
- adj: adjacency list storing edge indices per vertex
- n: number of vertices
- s: source vertex id
- t: sink vertex id
- level: BFS levels for Dinic
- ptr: current-edge iterators for DFS
- pi: node potentials to compute reduced costs rc = cost + pi[v] - pi[u]
- current_eps: current epsilon threshold; only edges with rc <= eps are considered
- level_adj: filtered adjacency by level and reduced-cost threshold
*/
struct Dinic{
    const long long flow_inf = (long long)4e18;
    vector<FlowEdge> edges;
    vector<vector<int>> adj;
    int n, m = 0;
    int s, t;
    vector<int> level, ptr;
    queue<int> q;

    vector<long long> pi;
    long double EPS_GROW = 2.0;
    long long EPS0 = 0;
    long long EPS_MAX = (long long)9e18;
    long long current_eps = 0;

    vector<vector<int>> level_adj;

    Dinic(int n, int s, int t) : n(n), s(s), t(t) {
        adj.resize(n);
        level.resize(n);
        ptr.resize(n);
        pi.assign(n, 0);
        level_adj.resize(n);
    }

    /*
    Insert a forward edge and its reverse edge into the residual network.
    - v: current node / tail node index.
    - u: neighbor / head node index.
    - cap: capacity of the edge.
    - cost: per-unit cost on the edge.
    */
    void add_edge(int v, int u, long long cap, long long cost = 0) {
        edges.emplace_back(v, u, cap, cost);
        edges.emplace_back(u, v, 0,  -cost);
        adj[v].push_back(m);
        adj[u].push_back(m + 1);
        m += 2;
    }

    inline long long rc(int id, int v, int u) const {
        __int128 val = (__int128)edges[id].cost + (__int128)pi[v] - (__int128)pi[u];
        if (val > (__int128)LLONG_MAX) return LLONG_MAX;
        if (val < (__int128)LLONG_MIN) return LLONG_MIN;
        return (long long)val;
    }
    
    /*
    Build level graph restricted to edges with rc <= current_eps to bias toward cheaper arcs.
    - MAXD: max BFS depth (level cutoff).
    */
    bool bfs(int MAXD) {
        // Reset levels and start BFS from source; only up to MAXD layers to bias short paths.
        fill(level.begin(), level.end(), -1);
        while(!q.empty()) q.pop();
        level[s] = 0; q.push(s);

        // Skip saturated edges and those with reduced cost rc > current_eps.
        while(!q.empty()){
        // BFS from all start-allowed vertices to compute shortest heads.
            int v = q.front(); q.pop();
            if (level[v] >= MAXD) continue;
            for (int id : adj[v]) {
                if (edges[id].cap == edges[id].flow) continue;
                int u = edges[id].u;
                if (rc(id, v, u) > current_eps) continue;
                if (level[u] != -1) continue;
                level[u] = level[v] + 1;
                q.push(u);
            }
        }
        if (level[t] == -1) return false;

        for (int v = 0; v < n; ++v) {
            // Build next-layer adjacency filtered by (level[u] == level[v]+1) and rc <= current_eps.
            level_adj[v].clear();
            if (level[v] == -1) continue;
            for (int id : adj[v]) {
                if (edges[id].cap == edges[id].flow) continue;
                int u = edges[id].u;
                if (level[u] != level[v] + 1) continue;
                if (rc(id, v, u) > current_eps) continue;
                level_adj[v].push_back(id);
            }
            sort(level_adj[v].begin(), level_adj[v].end(), [&](int a, int b){
                // Prefer cheaper reduced-cost arcs first; tie-break by larger residual capacity.
                int ua = edges[a].u, ub = edges[b].u;
                long long rca = rc(a, v, ua);
                long long rcb = rc(b, v, ub);
                if (rca != rcb) return rca < rcb;
                long long resa = edges[a].cap - edges[a].flow;
                long long resb = edges[b].cap - edges[b].flow;
                return resa > resb;
            });
        }
        return true;
    }
    
    /*
    Send blocking flow on the filtered level graph (Dinic DFS).
    - v: current node / tail node index.
    - pushed: flow amount we try to send.
    */
    long long dfs(int v, long long pushed) {
        // Base cases: nothing to push or sink reached.
        if (pushed == 0) return 0;
        if (v == t) return pushed;
        for (int &it = ptr[v]; it < (int)level_adj[v].size(); ++it) {
            // Iterate admissible edges in current level graph.
            int id = level_adj[v][it];
            int u = edges[id].u;
            long long avail = edges[id].cap - edges[id].flow;
            if (avail <= 0) continue;
            long long tr = dfs(u, (avail < pushed ? avail : pushed));
            if (tr == 0) continue;
            edges[id].flow += tr;
            edges[id ^ 1].flow -= tr;
            // Update residual capacities (forward/reverse edges).
            return tr;
        }
        return 0;
    }

    /*
    Find a max flow
    - MAXD: max BFS depth (level cutoff).
    */
    long long flow(int MAXD = (int)1e8) {
        // Epsilon-scheduled Dinic: tighten rc-threshold first; relax if we stall.
        long long f = 0;
        current_eps = EPS0;

        auto grow_eps = [&]()->bool {
            if (current_eps >= EPS_MAX) return false;
            if (current_eps <= 0) current_eps = 1;
            else {
                long double next = (long double)current_eps * EPS_GROW;
                if (next > (long double)EPS_MAX) current_eps = EPS_MAX;
                else current_eps = (long long)(next + 0.5L);
            }
            return true;
        };

        for (;;) {
            while (!bfs(MAXD)) {
            // Increase epsilon until a valid level graph is found or we max out.
                if (!grow_eps()) return f;
            }

            fill(ptr.begin(), ptr.end(), 0);
            // Push blocking flow over the filtered level graph.
            bool progressed = false;
            while (long long pushed = dfs(s, flow_inf)) {
                f += pushed;
                progressed = true;
            }

            if (!progressed) {
                // If no flow was pushed at this epsilon, relax threshold and retry.
                if (!grow_eps()) return f;
            }
        }
    }
};


/*
Builds an Euler tour over a multigraph. Assume that the Euler Tour exist
- G: multigraph adjacency
- E: collected trail triples
*/
struct Euler {
    vector<vector<ii>> G;
    vector<tuple<int,int,int>> E;

    Euler(const vector<vector<ii>>& grafo) : G(grafo) {}

    /*
    DFS style function to compute the Euler tour
    - u: node processing
    */
    void go(int u){
        // Consume edges; append in reverse order (we reverse once at the end).
        while(!G[u].empty()){
            auto [v, eid] = G[u].back(); G[u].pop_back();
            go(v);
            E.emplace_back(u, v, eid);
        }
    }

    /*
    Find Euler tour
    - u: node processing
    */
    vector<tuple<int,int,int>> get_path(int s){
        E.clear();
        go(s);
        reverse(E.begin(), E.end());
        return E;
    }
};


/*
A test action
- c: integer cost of executing this vertex
- loc: location constraint: 0=start, 1=mid, 2=end (3 is derived in Flow)
- pre: list of prerequisite vertex indices that must precede this vertex in a pipeline
*/
struct Vertex{
    ll c;
    int loc;
    vector<int> pre;
    Vertex(){}

    /*
    Read input.
    */
    void leer(){
        int Pl;
        cin>>c>>loc>>Pl;
        pre.resize(Pl);
        for(int j = 0; j<Pl; j++){
            cin>>pre[j];
        }
    }
};

/*
A flow (test segment): ordered vertex list with required multiplicity and derived endpoints/cost.
- c: aggregated cost of the flow (sum of vertex costs, respecting location rules)
- w: required number of times this flow must appear overall
- v: sequence of vertex ids in the flow
- s: start vertex id (first element of v)
- t: end vertex id (last element of v)
- loc: location type: 0=start-only, 2=end-only, 3=start+end, 1=free
*/
struct Flow{
    ll c = 0, w;
    vector<int> v;
    int s, t, loc;
    Flow(){}
    
    /*
    Read input and populate structure fields, deriving fast-to-use properties.
    - vertice: vertices vector.
    */
    void leer(vector<Vertex> &vertice){
        int L;
        cin>>w>>L;
        loc = 1;
        v.resize(L);
        for(int i = 0; i<L; i++){
            cin>>v[i];
            if(i == 0){
                s = v[i];
                if(vertice[v[i]].loc == 0){
                    loc = 0;
                    c += vertice[v[i]].c;
                }
                continue;
            }
            if(i == L-1){
                t = v[i];
                if(vertice[v[i]].loc == 2){
                    if(loc == 0){
                        loc = 3;
                    }
                    else{
                        loc = 2;
                    }
                }
            }
            c += vertice[v[i]].c;
        }
    }
};

/*
Select a minimal-cost subset of given pipelines preserving flow multiplicities and prerequisite coverage.
- pipeline: original pipelines (each is a list of flow ids)
- c: cost per original pipeline (sum of flows' costs)
- used: flag per original pipeline, 1 if kept in the subset
- w: covered multiplicity per flow id after selection
- req_pairs: universe of prerequisite pairs (p -> v) to cover
- pair2id: map (p,v) -> compact id used for coverage accounting
- fcnt: per-pipeline multiplicity of each flow id
- cov_pairs: per-pipeline list of (p->v) pair ids this pipeline covers
*/
struct Usub{
    int N;
    vector<vector<int>> pipeline;
    vector<long long> c;
    vector<char> used;
    vector<int> w;

    vector<pair<int,int>> req_pairs;
    unordered_map<uint64_t,int> pair2id;
    vector<vector<pair<int,int>>> fcnt;
    vector<vector<int>> cov_pairs;

    Usub(int n, int m, int k)
        : N(k),
          pipeline(k),
          c(k, 0),
          used(k, 0),
          w(m, 0) {}

    /*
    Read input and populate structure fields, deriving fast-to-use properties.
    -flow: flows vector.
    */
    void leer(const vector<Flow> &flow){
        // Read each pipeline and sum the costs of its flows.
        fill(used.begin(), used.end(), 0);
        for(int i = 0; i < N; i++){
            int L; cin >> L;
            pipeline[i].resize(L);
            long long cc = 0;
            for(int j = 0; j < L; j++){
                cin >> pipeline[i][j];
                int eid = pipeline[i][j];
                if (0 <= eid && eid < (int)flow.size()) cc += (long long)flow[eid].c;
            }
            c[i] = cc;
        }
    }

    vector<int> order_by_cost_desc() const {
        vector<int> ord(N); iota(ord.begin(), ord.end(), 0);
        sort(ord.begin(), ord.end(), [&](int a, int b){
            if (c[a] != c[b]) return c[a] > c[b];
            return a < b;
        });
        return ord;
    }

    /*
    Pack (a,b) into a 64-bit key for fast hashmap usage.
    */
    static inline uint64_t pack_pair_u64(uint32_t a, uint32_t b){
        return ( (uint64_t)a << 32 ) | (uint64_t)b;
    }

    /*
    Precompute per-pipeline multiplicities and prerequisite coverage (Usub/Unew flavor).
    -vertice: vertices vector.
    -flow: flows vector.
    */
    void pre(const vector<Vertex>& vertice,
             const vector<Flow>& flow)
    {
        // Build prerequisite pair universe (p->v) and source/target markings.
        const int n = (int)vertice.size();

        req_pairs.clear();
        pair2id.clear();
        req_pairs.reserve(1<<20);
        pair2id.reserve(1<<21);

        vector<char> isPreSrc(n, 0), isTarget(n, 0);
        for(int v=0; v<n; ++v){
            // Register each (p->v) and assign a compact id to the pair.
            if(!vertice[v].pre.empty()) isTarget[v] = 1;
            for(int p : vertice[v].pre){
                if(p>=0) isPreSrc[p] = 1;
                uint64_t key = pack_pair_u64((uint32_t)p,(uint32_t)v);
                if(!pair2id.count(key)){
                    int id = (int)req_pairs.size();
                    pair2id.emplace(key, id);
                    req_pairs.emplace_back(p, v);
                }
            }
        }
		
        fcnt.assign(N, {});
        cov_pairs.assign(N, {});
        // For each pipeline: (1) multiplicity of each flow via sorting+RLE; (2) which (p->v) pairs it covers.
        vector<int> seenMark(n, -1); int seenEpoch = 0;

        for(int i=0;i<N;i++){
            const auto &F = pipeline[i];
            if(!F.empty()){
                vector<int> tmp = F;
                // Count occurrences per flow id (sorted run-length encoding).
                sort(tmp.begin(), tmp.end());
                fcnt[i].reserve(tmp.size());
                for(size_t a=0; a<tmp.size(); ){
                    size_t b=a;
                    while(b<tmp.size() && tmp[b]==tmp[a]) ++b;
                    fcnt[i].push_back({tmp[a], (int)(b-a)});
                    a=b;
                }
            }

            vector<int> seq;
            if(!F.empty()){
            // Concatenate the vertex lists of flows; skip the first vertex of subsequent flows to avoid duplication at joints.
                seq.insert(seq.end(), flow[F[0]].v.begin(), flow[F[0]].v.end());
                for(size_t s=1; s<F.size(); ++s){
                    const auto &fv = flow[F[s]].v;
                    if(fv.size()>1) seq.insert(seq.end(), fv.begin()+1, fv.end());
                    else if(!fv.empty()) seq.push_back(fv[0]);
                }
            }

            vector<int> covered_ids; covered_ids.reserve(64);
            ++seenEpoch;
            for(int node : seq){
                // Mark seen prerequisites and record (p->v) coverage when all of P_v have been seen before v.
                if(isPreSrc[node]) seenMark[node] = seenEpoch;
                if(isTarget[node]){
                    for(int psrc : vertice[node].pre){
                        if(seenMark[psrc]==seenEpoch){
                            auto it = pair2id.find(pack_pair_u64((uint32_t)psrc,(uint32_t)node));
                            if(it!=pair2id.end()){
                                covered_ids.push_back(it->second);
                            }
                        }
                    }
                }
            }
            if(!covered_ids.empty()){
                sort(covered_ids.begin(), covered_ids.end());
                covered_ids.erase(unique(covered_ids.begin(), covered_ids.end()), covered_ids.end());
            }
            cov_pairs[i].swap(covered_ids);
        }
    }

    /*
    Drop expensive pipelines if quotas and prerequisite coverage remain satisfied.
    -flow: flows vector.
    */
    void greedy(const vector<Flow>& flow){
        // Try to drop costly pipelines while respecting quotas and prerequisite coverage.
        const int m = (int)flow.size();
        vector<char> keep(N, 1);

        vector<long long> coverFlow(m, 0);
        for(int i=0;i<N;i++) for(const auto &pr : fcnt[i]) coverFlow[pr.first] += pr.second;
        // coverFlow: total multiplicity per flow id across all pipelines.
        vector<long long> slack(m, 0);
        for(int j=0;j<m;j++) slack[j] = coverFlow[j] - flow[j].w;

        // slack tells how many occurrences we can remove per flow before violating its requirement.
        const int R = (int)req_pairs.size();
        vector<int> pairCover(R, 0);
        for(int i=0;i<N;i++) for(int id : cov_pairs[i]) pairCover[id]++;

        // pairCover: how many pipelines currently cover each (p->v) pair.
        auto removable = [&](int i)->bool{
            // Removable if all covered pairs remain covered and slack is sufficient for all of its flows.
            for(int id : cov_pairs[i]) if(pairCover[id] <= 1) return false;
            for(const auto &pr : fcnt[i]){
                int fid = pr.first, cnt = pr.second;
                if(slack[fid] < cnt) return false;
            }
            return true;
        };

        auto ord = order_by_cost_desc();

        // Process pipelines from most to least expensive.
        for(int i : ord){
            // Pass 1: greedy removals in cost order.
            if(!keep[i]) continue;
            if(removable(i)){
                keep[i]=0;
                for(const auto &pr : fcnt[i]){
                    int fid=pr.first, cnt=pr.second;
                    slack[fid] -= cnt;
                    coverFlow[fid] -= cnt;
                }
                for(int id : cov_pairs[i]) pairCover[id]--;
            }
        }
        for(int i=0;i<N;i++){
            if(!keep[i]) continue;
        // Pass 2: another sweep may become possible after adjustments.
            if(removable(i)){
                keep[i]=0;
                for(const auto &pr : fcnt[i]){
                    int fid=pr.first, cnt=pr.second;
                    slack[fid] -= cnt;
                    coverFlow[fid] -= cnt;
                }
                for(int id : cov_pairs[i]) pairCover[id]--;
            }
        }

        fill(used.begin(), used.end(), 0);
        fill(w.begin(), w.end(), 0);
        for(int i=0;i<N;i++){
            if(!keep[i]) continue;
            used[i] = 1;
            for(auto &pr : fcnt[i]) w[pr.first] += pr.second;
        }
    }

    /*
    Execute the full algorithm for this subproblem.
    -vertice: vertices vector.
    -flow: flows vector.
    */
    void run(vector<Vertex> &vertice, vector<Flow> &flow){
        pre(vertice, flow);
        greedy(flow);
    }

    /*
    Print the component’s required output format.
    */
    void imprimir(){
        int L = 0; for(char u : used) L += (int)u;
        cout << L << "\n";
        for(int i = 0; i < N; i++) if (used[i]) cout << i << " ";
        cout << "\n";
    }
};

/*
Construct new pipelines, from via degree balancing, cost-guided max-flow and Euler trails.
- pipeline: constructed pipelines (each is a list of flow ids)
- G: auxiliary multigraph (u -> v, with edgeId for flow or -1 as separator)
- degree: balance per vertex to ensure concatenability (in-degree vs out-degree)
- w: multiplicity achieved per flow id after construction
*/
struct Unew{
    int N, M;
    vector<vector<int> > pipeline;
    vector<vector<ii> > G;
    vector<int> degree;
    vector<int> w;

    Unew(int n, int m) : N(n), M(m), G(n+1), degree(n+1) {}

    /*
    Recompute multiplicity per flow id from current pipelines.
    -flow: flows vector.
    */
    void computew(vector<Flow> &flow){
		w.clear();
		w.resize(M, 0);
		for(auto &u: pipeline){
			for(auto &&uu: u){
				w[uu]++;
			}
		}
	}

    /*
    Precompute per-pipeline multiplicities and prerequisite coverage (Usub/Unew flavor).
    -U: given pipelines in input (Ufull).
    -flow: flows vector.
    -vertice: vertices vector.
    */
    void pre(const Usub& U, const vector<Flow>& flow, const vector<Vertex>& vertice){
        // Seed new set with pipelines that add prerequisite coverage relative to current state.
        const int n = (int)vertice.size();
        const int K = U.N;

        vector<char> isTarget(n, 0), isPreSrc(n, 0), covered(n, 0);
        for(int v = 0; v<n; v++){
            if(vertice[v].pre.empty()) covered[v] = 1;
            else isTarget[v] = 1;
            for(int p : vertice[v].pre) if(p>=0 && p<n) isPreSrc[p] = 1;
        }

        vector<int> seenMark(n, -1), gotMark(n, -1);
        int epoch = 0;

        auto cover_from_pipeline = [&](const vector<int>& flist){
            // Simulate concatenated vertices of flist; mark targets whose prerequisites already appeared.
            if(flist.empty()) return;
            vector<int> seq;
            seq.insert(seq.end(), flow[flist[0]].v.begin(), flow[flist[0]].v.end());
            for(size_t s=1; s<flist.size(); ++s){
                const auto &vv = flow[flist[s]].v;
                if (vv.size() > 1) seq.insert(seq.end(), vv.begin()+1, vv.end());
                else if(!vv.empty()) seq.push_back(vv[0]);
            }
            ++epoch;
            for(int node : seq){
                if(isTarget[node] && gotMark[node] != epoch){
                    bool ok = true;
                    const auto &P = vertice[node].pre;
                    for(int p : P){
                        if (p < 0 || p >= n) { ok = false; break; }
                        if (seenMark[p] != epoch) { ok = false; break; }
                    }
                    if(ok){
                        gotMark[node] = epoch;
                        covered[node] = 1;
                    }
                }
                if(isPreSrc[node]) seenMark[node] = epoch;
            }
        };

        for(const auto &flist : pipeline) cover_from_pipeline(flist);

        vector<int> ord(K); iota(ord.begin(), ord.end(), 0);
        stable_sort(ord.begin(), ord.end(), [&](int a, int b){        // Favor lower-cost pipelines when seeding (stable among equals).

            if (U.c[a] != U.c[b]) return U.c[a] < U.c[b];
            return a < b;
        });

        for(int pid : ord){
            const auto &flist = U.pipeline[pid];
            if(flist.empty()) continue;

            vector<int> seq;
            seq.insert(seq.end(), flow[flist[0]].v.begin(), flow[flist[0]].v.end());
            for(size_t s=1; s<flist.size(); ++s){
                const auto &vv = flow[flist[s]].v;
                if (vv.size() > 1) seq.insert(seq.end(), vv.begin()+1, vv.end());
                else if(!vv.empty()) seq.push_back(vv[0]);
            }

            ++epoch;
            int gain = 0;
            vector<int> newly; newly.reserve(16);

            for(int node : seq){
                if(isTarget[node] && gotMark[node] != epoch){
                    bool ok = true;
                    const auto &P = vertice[node].pre;
                    for(int p : P){
                        if (p < 0 || p >= n) { ok = false; break; }
                        if (seenMark[p] != epoch) { ok = false; break; }
                    }
                    if(ok){
                        gotMark[node] = epoch;
                        if(!covered[node]) { gain++; newly.push_back(node); }
                    }
                }
                if(isPreSrc[node]) seenMark[node] = epoch;
            }

            if(gain > 0){
                // Keep this pipeline; it adds at least one newly covered target.
                pipeline.push_back(flist);
                for(int v : newly) covered[v] = 1;
            }
        }
    }

    /*
    Expand remaining quotas into a multigraph with separators by location rules.
    -flow: flows vector.
    */
    void makegraph(vector<Flow> &flow){
        // Convert remaining flow quotas into degree balances and auxiliary edges (with -1 separators for location rules).
        w.clear();
        w.resize(flow.size(), 0);
        for(auto u: pipeline){
            for(auto uu: u){
                w[uu]++;
            }
        }
        for(int i = 0; i<M; i++){
            for(int j = 0; j<flow[i].w - w[i]; j++){                // For each missing occurrence of flow i, add s->t; also add N<->s/t edges if start/end constraints apply.

                if(flow[i].loc == 0 || flow[i].loc == 3){
                    degree[N]++;
                    degree[flow[i].s]--;
                    G[N].push_back({flow[i].s, -1});
                }
                if(flow[i].loc == 2 || flow[i].loc == 3){
                    degree[flow[i].t]++;
                    degree[N]--;
                    G[flow[i].t].push_back({N, -1});
                }
                degree[flow[i].s]++;
                degree[flow[i].t]--;
                G[flow[i].s].push_back({flow[i].t, i});
            }
        }
    }

    /*
    Use cost-guided Dinic to connect surpluses/deficits in degrees and record chosen transitions.
    -flow: flows vector.
    */
    void completegraph(vector<Flow> &flow){
        // Satisfy degree imbalances via super source/sink; connect s->t arcs by increasing cost using epsilon-filtered Dinic.
        Dinic dinic(N+3, N+1, N+2);
        int needed = 0;
        for(int i = 0; i<=N; i++){
            // degree<0: need inflow from super-source (N+1); degree>0: send to super-sink (N+2).
            if(degree[i] < 0){
                dinic.add_edge(N+1, i, -degree[i]);
            }
            if(degree[i] > 0){
                needed += degree[i];
                dinic.add_edge(i, N+2, degree[i]);
            }
        }

        map<ii, ii> arista;

        for(int i = 0; i<M; i++){
            Flow u = flow[i];            // Allow jumps from sentinel N to starts, and from ends to N, so we can open/close pipelines.

            if(u.loc == 0 || u.loc == 3){
                dinic.add_edge(N, u.s, 1e9);
                arista[{N, u.s}] = {0, -1};
            }
            if(u.loc == 2 || u.loc == 3){
                dinic.add_edge(u.t, N, 1e9);
                arista[{u.t, N}] = {0, -1};
            }
        }

        vector<vector<ll> > orden(M, vector<ll>(5));

        for(int i = 0; i<M; i++){
            orden[i][0] = flow[i].c;
            orden[i][1] = flow[i].s;
            orden[i][2] = flow[i].t;
            orden[i][3] = i;
            orden[i][4] = flow[i].w;
        }

        sort(orden.begin(), orden.end());
        // Sort candidate arcs by flow cost; we will add them in three stages (cheap→expensive).
		int limit = M/3;
		int limit2 = 2*M/3;
		dinic.EPS_GROW = 2.0;
    	dinic.EPS0 = ll(1e18);
    	dinic.EPS_MAX = -1;

        for(int i = 0; i<limit; i++){
            // Stage A: cheapest third, very tight epsilon.
            int w = orden[i][4];
            int s = orden[i][1];
            int t = orden[i][2];
            int idx = orden[i][3];
            if(arista.count({s, t}) == 0){
                arista[{s, t}] = {w, idx};
            }
            dinic.add_edge(s, t, 1e9, orden[i][0]);
			dinic.EPS0 = min(dinic.EPS0, orden[i][0]);
			dinic.EPS_MAX = max(dinic.EPS_MAX, orden[i][0]);
        }
		dinic.EPS_GROW = sqrt(sqrt(dinic.EPS_MAX/dinic.EPS0)) + 2;
        dinic.flow(20);
        // Short run to exploit cheapest opportunities before expanding the candidate set.
		for(int i = limit; i<limit2; i++){
            // Stage B: mid third, slightly relaxed epsilon.
            int w = orden[i][4];
            int s = orden[i][1];
            int t = orden[i][2];
            int idx = orden[i][3];
            if(arista.count({s, t}) == 0){
                arista[{s, t}] = {w, idx};
            }
            dinic.add_edge(s, t, 1e9, orden[i][0]);
			dinic.EPS0 = min(dinic.EPS0, orden[i][0]);
			dinic.EPS_MAX = max(dinic.EPS_MAX, orden[i][0]);
        }
		dinic.EPS_GROW = sqrt(sqrt(sqrt(dinic.EPS_MAX/dinic.EPS0))) + 2;
		dinic.flow(10);
		for(int i = limit2; i<M; i++){
            int w = orden[i][4];
            int s = orden[i][1];
            int t = orden[i][2];
            int idx = orden[i][3];
            if(arista.count({s, t}) == 0){
                arista[{s, t}] = {w, idx};
            }
            dinic.add_edge(s, t, 1e9, orden[i][0]);
			dinic.EPS0 = min(dinic.EPS0, orden[i][0]);
			dinic.EPS_MAX = max(dinic.EPS_MAX, orden[i][0]);
        }
		dinic.EPS_GROW = sqrt(sqrt(dinic.EPS_MAX/dinic.EPS0)) + 2;
		dinic.flow();
        
        for(int id = 0; id<int(dinic.edges.size()); id++){
        // Materialize positive-flow arcs into multigraph G (repeat per unit of flow).
            auto &E = dinic.edges[id];
            if(E.flow <= 0) continue;

            int v = E.v, u = E.u;
            if (v == N+1 || u == N+2) continue;

            for(ll f = 0; f<E.flow; f++){
                degree[v]++; degree[u]--;
                G[v].push_back({u, arista[{v, u}].second});
            }
        }
    }

    /*
    Extract pipelines by Euler trail; -1 edges split sequences into separate pipelines.
    -flow: flows vector.
    */
    void buildpipelines(vector<Flow> &flow){

        // Shuffle adjacency to reduce bias; then derive Euler trail from sentinel N and split by -1 separators.
        for(auto &u: G){
            random_shuffle(u.begin(), u.end());
        }

        Euler euler(G);
        vector<tuple<int,int,int>> P = euler.get_path(N);

        vector<int> cur;
        for(auto u: P){
            auto [from, to, eid] = u;
            if(eid == -1){
            // Separator: end current pipeline and start a new one.
                if(!cur.empty()){
                    pipeline.push_back(cur);
                    cur.clear();
                }
            }
            else{
                cur.push_back(eid);
            // Append flow id to current pipeline.
            }
        }
        if (!cur.empty()) pipeline.push_back(cur);
        
    }

    struct Edge{
        int to;
        ll dis;
        int idx;
        Edge(int to, ll dis, int idx) : to(to), dis(dis), idx(idx) {}
    };

    /*
    If a pipeline is too long, split using short forward/backward chains around a boundary vertex.
    -flow: flows vector.
    -MAX_SEG: maximum long parameter.
    */
    void rebalance(vector<Flow> &flow, int MAX_SEG = 1000){
        // Split pipelines exceeding MAX_SEG by bridging to nearest end and start via cheap chains.
        struct E { int to, eid; };
        const int V = N;
        const int INF = 1e9;

        vector<vector<E>> g(V), rg(V);
		g.reserve(V); rg.reserve(V);

		for (int eid = 0; eid < (int)flow.size(); ++eid) {
            // Build forward (g) and reverse (rg) adjacency lists from flow endpoints.
			int a = flow[eid].s, b = flow[eid].t;
			if (0 <= a && a < V && 0 <= b && b < V) {
				g[a].push_back({b, eid});
				rg[b].push_back({a, eid});
			}
		}

		auto by_cost = [&](const E& A, const E& B){
			if (flow[A.eid].c != flow[B.eid].c) return flow[A.eid].c < flow[B.eid].c;
			return A.eid < B.eid;
		};

		for (int u = 0; u < V; ++u) {
			sort(g[u].begin(),  g[u].end(),  by_cost);
        // Prefer cheaper connector edges first.
			sort(rg[u].begin(), rg[u].end(), by_cost);
		}

        vector<int> distF(V, INF);
        vector<pair<int,int>> prevF(V, {-1, -1});
        deque<int> q;
        for (int eid = 0; eid < (int)flow.size(); ++eid) {
            if (flow[eid].loc == 0 || flow[eid].loc == 3) {
                int s = flow[eid].s;
                if (0 <= s && s < V && distF[s] > 0) {
                    distF[s] = 0;
                    prevF[s] = {-1, -1};
                    q.push_back(s);
                }
            }
        }
        while (!q.empty()) {
            int v = q.front(); q.pop_front();
            int dv = distF[v];
            for (const auto &e : g[v]) {
                if (distF[e.to] > dv + 1) {
                    distF[e.to] = dv + 1;
                    prevF[e.to] = {v, e.eid};
                    q.push_back(e.to);
                }
            }
        }

        vector<int> distR(V, INF);
        vector<pair<int,int>> nextR(V, {-1, -1});
        deque<int> rq;
        for (int eid = 0; eid < (int)flow.size(); ++eid) {
            if (flow[eid].loc == 2 || flow[eid].loc == 3) {
                int t = flow[eid].t;
                if (0 <= t && t < V && distR[t] > 0) {
                    distR[t] = 0;
                    nextR[t] = {-1, -1};
                    rq.push_back(t);
                }
            }
        }
        while (!rq.empty()) {
        // Reverse BFS from all end-allowed vertices to compute shortest tails.
            int v = rq.front(); rq.pop_front();
            int dv = distR[v];
            for (const auto &e : rg[v]) {
                if (distR[e.to] > dv + 1) {
                    distR[e.to] = dv + 1;
                    nextR[e.to] = {v, e.eid};
                    rq.push_back(e.to);
                }
            }
        }

        auto boundary_vertex = [&](const vector<int>& P, int j)->int{
            return flow[P[j]].t;
        };

        for (int i = 0; i < (int)pipeline.size(); ++i) {
            auto &P = pipeline[i];
            while ((int)P.size() > MAX_SEG) {
                // Try to cut near MAX_SEG; ensure the suffix can be completed to an end and the prefix from a start.
                const int L = (int)P.size();
                int limit = min(MAX_SEG - 1, L - 1);
                int cut = -1;

                for (int j = limit; j >= 0; --j) {
                    // Scan candidate cut positions from right to left.
                    int vcut = boundary_vertex(P, j);
                    if (vcut < 0 || vcut >= V) continue;
                    if (distR[vcut] == INF || distF[vcut] == INF) continue;
                    if ((j + 1) + distR[vcut] <= MAX_SEG) {
                        cut = j; break;
                    }
                }
                if (cut == -1) break;

                const int vcut = boundary_vertex(P, cut);

                vector<int> suffix;
                suffix.reserve(L - (cut + 1));
                suffix.insert(suffix.end(), P.begin() + (cut + 1), P.end());

                P.resize(cut + 1);
                {
                    int x = vcut;
                    // Append short tail from vcut to a valid end using reverse parents.
                    while (nextR[x].first != -1) {
                        P.push_back(nextR[x].second);
                        x = nextR[x].first;
                    }
                }

                vector<int> pre;
                {
                    int y = vcut;
                    // Prepend short head from a valid start to vcut using forward parents.
                    while (prevF[y].first != -1) {
                        pre.push_back(prevF[y].second);
                        y = prevF[y].first;
                    }
                    reverse(pre.begin(), pre.end());
                }

                pipeline.emplace_back();
                auto &Pnew = pipeline.back();
                Pnew.swap(suffix);
                Pnew.insert(Pnew.begin(), pre.begin(), pre.end());
            }
        }
    }

    /*
    Optional pass to remove redundant pipelines while keeping quotas/coverage.
    -flow: flows vector.
    -vertice: vertices vector.
    */
	void prune(const vector<Flow>& flow, const vector<Vertex>& vertice){
        // Recompute coverage and multiplicities on constructed pipelines; drop those that are safe to remove.
		const int K = (int)pipeline.size();
		if (K == 0) { w.assign(flow.size(), 0); return; }

		const int n = (int)vertice.size();
		const int m = (int)flow.size();

		auto pack_u64 = [](uint32_t a, uint32_t b)->uint64_t {
			return ( (uint64_t)a << 32 ) | (uint64_t)b;
		};
		unordered_map<uint64_t,int> pair2id; pair2id.reserve(1<<21);
		vector<pair<int,int>> req_pairs;     req_pairs.reserve(1<<20);
		vector<char> isPreSrc(n, 0), isTarget(n, 0);
		for(int v=0; v<n; ++v){
			if(!vertice[v].pre.empty()) isTarget[v] = 1;
			for(int p : vertice[v].pre){
				if(p>=0 && p<n) isPreSrc[p] = 1;
				uint64_t key = pack_u64((uint32_t)max(0,p), (uint32_t)v);
				if(!pair2id.count(key)){
					int id = (int)req_pairs.size();
					pair2id.emplace(key, id);
					req_pairs.emplace_back(p, v);
				}
			}
		}
		const int R = (int)req_pairs.size();

		vector<vector<pair<int,int>>> fcnt(K);
		vector<vector<int>> cov(K);
		vector<long long> pcost(K, 0);

		vector<int> seen(n, -1); int epoch = 0;

		for(int i=0;i<K;i++){
			const auto &F = pipeline[i];
			if(!F.empty()){
				pcost[i] = 0;
				vector<int> tmp = F;
				sort(tmp.begin(), tmp.end());
				fcnt[i].reserve(tmp.size());
				for(size_t a=0; a<tmp.size(); ){
					size_t b=a;
					while(b<tmp.size() && tmp[b]==tmp[a]) ++b;
					int fid = tmp[a], cnt = (int)(b-a);
					fcnt[i].push_back({fid, cnt});
					if (0 <= fid && fid < m) pcost[i] += (long long)flow[fid].c * cnt;
					a=b;
				}
			}

			vector<int> seq;
			if(!F.empty()){
				seq.insert(seq.end(), flow[F[0]].v.begin(), flow[F[0]].v.end());
				for(size_t s=1; s<F.size(); ++s){
					const auto &vv = flow[F[s]].v;
					if(vv.size()>1) seq.insert(seq.end(), vv.begin()+1, vv.end());
					else if(!vv.empty()) seq.push_back(vv[0]);
				}
			}

			vector<int> covered_ids; covered_ids.reserve(32);
			++epoch;
			for(int node : seq){
				if(isTarget[node]){
					for(int psrc : vertice[node].pre){
						if(psrc>=0 && psrc<n && seen[psrc]==epoch){
							auto it = pair2id.find(pack_u64((uint32_t)psrc,(uint32_t)node));
							if(it!=pair2id.end()) covered_ids.push_back(it->second);
						}
					}
				}
				if(isPreSrc[node]) seen[node] = epoch;
			}
			if(!covered_ids.empty()){
				sort(covered_ids.begin(), covered_ids.end());
				covered_ids.erase(unique(covered_ids.begin(), covered_ids.end()), covered_ids.end());
			}
			cov[i].swap(covered_ids);
		}

		vector<long long> coverFlow(m, 0);
		for(int i=0;i<K;i++) for(const auto &pr : fcnt[i]){
			if (0 <= pr.first && pr.first < m) coverFlow[pr.first] += pr.second;
		}
		vector<int> pairCover(R, 0);
		for(int i=0;i<K;i++) for(int id : cov[i]) pairCover[id]++;

		auto canDrop = [&](int i)->bool{
			for(int id : cov[i]) if(pairCover[id] <= 1) return false;
			for(const auto &pr : fcnt[i]){
				int fid = pr.first, cnt = pr.second;
				if(fid < 0 || fid >= m) continue;
				if(coverFlow[fid] - cnt < flow[fid].w) return false;
			}
			return true;
		};

		vector<int> ord(K); iota(ord.begin(), ord.end(), 0);
		sort(ord.begin(), ord.end(), [&](int a, int b){
			if (pcost[a] != pcost[b]) return pcost[a] > pcost[b];
			return (int)pipeline[a].size() > (int)pipeline[b].size();
		});

		vector<char> keep(K, 1);
		for(int i : ord){
			if(!keep[i]) continue;
			if(canDrop(i)){
				keep[i] = 0;
				for(int id : cov[i]) pairCover[id]--;
				for(const auto &pr : fcnt[i]){
					int fid = pr.first, cnt = pr.second;
					if(0 <= fid && fid < m) coverFlow[fid] -= cnt;
				}
			}
		}

		for(int i=0;i<K;i++){
			if(!keep[i]) continue;
			if(canDrop(i)){
				keep[i] = 0;
				for(int id : cov[i]) pairCover[id]--;
				for(const auto &pr : fcnt[i]){
					int fid = pr.first, cnt = pr.second;
					if(0 <= fid && fid < m) coverFlow[fid] -= cnt;
				}
			}
		}

		vector<vector<int>> np; np.reserve(K);
		for(int i=0;i<K;i++) if(keep[i]) np.push_back(std::move(pipeline[i]));
		pipeline.swap(np);

		w.assign(m, 0);
		for(int i=0;i<(int)pipeline.size(); i++){
			for(int fid : pipeline[i]){
				if(0 <= fid && fid < m) w[fid] += 1;
			}
		}
	}

    /*
    Verify that every flow multiplicity meets its required quota.
    -flow: flows vector.
    */
	bool check(vector<Flow> &flow){
		computew(flow);
		for(int i = 0; i<M; i++){
			if(w[i] < flow[i].w) return false;
		}
		return true;
	}

    /*
    Execute the full algorithm for this component.
    -U: given pipelines Ufull.
    -flow: flows vector.
    -vertice: vertices vector.
    */
    void run(Usub U, vector<Flow> &flow, vector<Vertex> &vertice){
        pre(U, flow, vertice);
		makegraph(flow);
        completegraph(flow);
        buildpipelines(flow);
        rebalance(flow);
		if(check(flow) == false){
			pipeline.clear();
			G.clear();
			degree.clear();
			G.resize(N+1);
			degree.resize(N+1);
			makegraph(flow);
			completegraph(flow);
			buildpipelines(flow);
			rebalance(flow);
			pre(U, flow, vertice);
		}

    }

    /*
    Print the component’s required output format.
    */
    void imprimir(){
        cout<<int(pipeline.size())<<"\n";
        for(auto &u: pipeline){
            cout<<int(u.size())<<" ";
            for(auto &uu: u){
                cout<<uu<<" ";
            }
            cout<<"\n";
        }
    }
};

int main(){

    // Driver: read input, solve U_sub, print; then build U_new and print.
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n, m, k;
    cin>>n>>m>>k;

    vector<Flow> flow(m);
    vector<Vertex> vertice(n);
    Usub U(n, m, k);

    vertice.resize(n);
    for(int i = 0; i<n; i++){
        vertice[i].leer();
    }

    for(int i = 0; i<m; i++){
        flow[i].leer(vertice);
    }

    U.leer(flow);
    U.run(vertice, flow);
    U.imprimir();

    Unew V(n, m);

    V.run(U, flow, vertice);

    V.imprimir();

    return 0;
}