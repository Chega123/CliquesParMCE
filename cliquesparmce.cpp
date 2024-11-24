#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <tbb/tbb.h>
#include <tbb/concurrent_vector.h>
#include <tbb/concurrent_unordered_set.h>

using Graph = std::unordered_map<int, std::unordered_set<int>>;
using Set = std::unordered_set<int>;
using namespace std;

struct SetHash {
    size_t operator()(const Set &s) const {
        size_t hash = 0;
        for (int v : s) {
            hash ^= std::hash<int>{}(v) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        }
        return hash;
    }
};

Graph loadGraph(const std::string &filename, int &numNodes, int &numEdges) {
    Graph G;
    std::ifstream infile(filename);
    int u, v;

    while (infile >> u >> v) {
        if (G[u].insert(v).second) {
            numEdges++;
        }
        G[v].insert(u);
    }

    numNodes = G.size();
    return G;
}

int graph_rank(const Graph &G, int v) {
    return G.at(v).size();
}

int ParPivot(const Graph &G, const Set &cand) {
    int pivot = -1;
    int maxDegree = -1;
    for (int v : cand) {
        int degree = G.at(v).size();
        if (degree > maxDegree) {
            maxDegree = degree;
            pivot = v;
        }
    }
    return pivot;
}

void ParTTT(const Graph &G, Set K, Set cand, Set fini, tbb::concurrent_unordered_set<Set, SetHash> &maximalCliques) {
    if (cand.empty() && fini.empty()) {
        std::vector<int> sortedK(K.begin(), K.end());
        std::sort(sortedK.begin(), sortedK.end());
        Set sortedSet(sortedK.begin(), sortedK.end());
        maximalCliques.insert(sortedSet);
        return;
    }

    int pivot = ParPivot(G, cand);
    Set ext;

    for (int v : cand) {
        if (G.at(pivot).find(v) == G.at(pivot).end()) {
            ext.insert(v);
        }
    }

    tbb::parallel_for(size_t(0), ext.size(), [&](size_t i) {
        auto it = std::next(ext.begin(), i);
        int q = *it;
        
        Set Kq = K;
        Kq.insert(q);

        Set candq, finiq;
        for (int neighbor : G.at(q)) {
            if (cand.find(neighbor) != cand.end()) {
                candq.insert(neighbor);
            }
            if (fini.find(neighbor) != fini.end()) {
                finiq.insert(neighbor);
            }
        }

        if (candq.empty() && finiq.empty()) {
            std::vector<int> sortedKq(Kq.begin(), Kq.end());
            std::sort(sortedKq.begin(), sortedKq.end());
            Set sortedSet(sortedKq.begin(), sortedKq.end());
            maximalCliques.insert(sortedSet);
        } else {
            ParTTT(G, Kq, candq, finiq, maximalCliques);
        }
    });
}

void ParMCE(const Graph &G, tbb::concurrent_unordered_set<Set, SetHash> &maximalCliques) {
    tbb::parallel_for_each(G.begin(), G.end(), [&](const std::pair<const int, std::unordered_set<int>> &entry) {
    //tbb::parallel_for_each(G.begin(), G.end(), [&](const auto &entry) {
        int v = entry.first;
        
        Set neighbors = G.at(v);
        neighbors.insert(v);
        Graph Gv;
        for (int u : neighbors) {
            for (int w : G.at(u)) {
                if (neighbors.find(w) != neighbors.end()) {
                    Gv[u].insert(w);
                }
            }
        }

        Set K = {v};
        Set cand, fini;

        for (int w : G.at(v)) {
            if (graph_rank(G, w) > graph_rank(G, v)) {
                cand.insert(w);
            } else {
                fini.insert(w);
            }
        }

        ParTTT(Gv, K, cand, fini, maximalCliques);
    });
}

int main() {
    int numNodes = 0, numEdges = 0;
    string dataset = "CA-GrQc.txt";
    Graph G = loadGraph(dataset, numNodes, numEdges);
    int num_threads = 16;

    cout << "Dataset: " << dataset << endl;
    cout << "Number of nodes: " << numNodes << endl;
    cout << "Number of edges: " << numEdges << endl;

    tbb::global_control control(tbb::global_control::max_allowed_parallelism, num_threads);

    tbb::concurrent_unordered_set<Set, SetHash> maximalCliques;
    
    tbb::tick_count start = tbb::tick_count::now();
    ParMCE(G, maximalCliques);
    tbb::tick_count end = tbb::tick_count::now();
    double elapsedTime = (end - start).seconds();

    cout << "Time taken: " << elapsedTime << " seconds" << endl;
    cout << "Total maximal cliques found: " << maximalCliques.size() << endl;

    return 0;
}
