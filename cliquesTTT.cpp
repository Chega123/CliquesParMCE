#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>


using Graph = unordered_map<int, unordered_set<int>>;
using Set = unordered_set<int>;
using namespace std;

Graph loadGraph(const string &filename, int &numNodes, int &numEdges) {
    Graph G;
    ifstream infile(filename);
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

// Función para seleccionar el pivote con el mayor número de vecinos en cand ∪ fini
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


void TTT(const Graph &G, Set K, Set cand, Set fini, vector<Set> &maximalCliques) {
    if (cand.empty() && fini.empty()) {
        maximalCliques.push_back(K);  // agrega k a la lista de cliques maximales
        return;
    }
    int pivot = Pivot(G, cand, fini);    // seleccion del pivote
    Set ext;    // creacion del conjunto ext 
    for (int v : cand) {
        if (G.at(pivot).find(v) == G.at(pivot).end()) {
            ext.insert(v);
        }
    }

    // procesamos cada vértice ext
    for (int q : ext) {
        Set Kq = K;
        Kq.insert(q);

        // construir candq y finiq
        Set candq, finiq;
        for (int neighbor : G.at(q)) {
            if (cand.find(neighbor) != cand.end()) {
                candq.insert(neighbor);
            }
            if (fini.find(neighbor) != fini.end()) {
                finiq.insert(neighbor);
            }
        }
        // actualizar cand y fini para eliminar q
        cand.erase(q);
        fini.insert(q);
        TTT(G, Kq, candq, finiq, maximalCliques);
    }
}


int main() {
    int numNodes = 0, numEdges = 0;
    string dataset = "/home/chega/cliquesparalela/dataset/CA-AstroPh.txt";
    Graph G = loadGraph(dataset, numNodes, numEdges);
    cout << "Dataset: " << dataset << endl;

    cout << "Number of nodes: " << numNodes << endl;
    cout << "Number of edges: " << numEdges << endl;

    Set K, cand, fini;
    for (const auto &[node, _] : G) {
        cand.insert(node);
    }

    vector<Set> maximalCliques;
    clock_t start = clock();
    TTT(G, K, cand, fini, maximalCliques);
    clock_t end = clock();
    double elapsedTime = double(end - start) / CLOCKS_PER_SEC;
    cout << "Time taken: " << elapsedTime << " seconds" << endl;
    cout << "Total maximal cliques found: " << maximalCliques.size() << endl;

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
