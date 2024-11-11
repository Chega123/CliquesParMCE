#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>

using Graph = std::unordered_map<int, std::unordered_set<int>>;
using Set = std::unordered_set<int>;
using namespace std;

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

// para seleccionar el pivote con el mayor número de vecinos en cand ∪ fini
int Pivot(const Graph &G, const Set &cand, const Set &fini) {
    int pivot = -1;
    int maxDegree = -1;
    for (int v : cand) {
        int degree = G.at(v).size();
        if (degree > maxDegree) {
            maxDegree = degree;
            pivot = v;
        }
    }
    for (int v : fini) {
        int degree = G.at(v).size();
        if (degree > maxDegree) {
            maxDegree = degree;
            pivot = v;
        }
    }
    return pivot;
}

// TTT secuencial
void TTT(const Graph &G, Set K, Set cand, Set fini, vector<Set> &maximalCliques) {
    if (cand.empty() && fini.empty()) {
        maximalCliques.push_back(K);  // agrega K a la lista de cliques maximales
        return;
    }

    // escojo pivote
    int pivot = Pivot(G, cand, fini);

    // creo conjunto ext = cand - ΓG(pivot)
    Set ext;
    for (int v : cand) {
        if (G.at(pivot).find(v) == G.at(pivot).end()) {
            ext.insert(v);
        }
    }

    // procesa cada vertice q en ext
    for (int q : ext) {
        Set Kq = K;
        Kq.insert(q);

        // candq y finiq
        Set candq, finiq;
        for (int neighbor : G.at(q)) {
            if (cand.find(neighbor) != cand.end()) {
                candq.insert(neighbor);
            }
            if (fini.find(neighbor) != fini.end()) {
                finiq.insert(neighbor);
            }
        }

        // actualiza cand y fini para eliminar q
        cand.erase(q);
        fini.insert(q);

        // recursivo
        TTT(G, Kq, candq, finiq, maximalCliques);
    }
}


int main() {
    int numNodes = 0, numEdges = 0;
    string dataset = "/home/chega/cliquesparalela/dataset/CA-AstroPh.txt";
    Graph G = loadGraph(dataset, numNodes, numEdges);
    //Graph G = loadGraph("com-dblp.ungraph.txt", numNodes, numEdges);
    //Graph G = loadGraph("dataset_test.txt", numNodes, numEdges);
    std::cout << "Dataset: " << dataset << std::endl;

    std::cout << "Number of nodes: " << numNodes << std::endl;
    std::cout << "Number of edges: " << numEdges << std::endl;

    Set K, cand, fini;
    for (const auto &[node, _] : G) {
        cand.insert(node);
    }

    vector<Set> maximalCliques;

    // TTT secuencial
    clock_t start = clock();
    TTT(G, K, cand, fini, maximalCliques);
    clock_t end = clock();
    double elapsedTime = double(end - start) / CLOCKS_PER_SEC;
    std::cout << "Time taken: " << elapsedTime << " seconds" << std::endl;
    std::cout << "Total maximal cliques found: " << maximalCliques.size() << std::endl;

    /*
    cout << "Maximal cliques:" << endl;
    for (const auto &clique : maximalCliques) {
        for (int v : clique) {
            cout << v << " ";
        }
        cout << endl;
    }
    */
   
    return 0;
}
