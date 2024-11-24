#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <tbb/tbb.h>
#include <tbb/concurrent_unordered_map.h>
#include <tbb/parallel_for_each.h>
#include <tbb/concurrent_vector.h>
#include <tbb/parallel_for.h>
#include <tbb/mutex.h>

// Type definitions for ease
using Graph = std::unordered_map<int, std::unordered_set<int>>;
using Set = std::unordered_set<int>;
using namespace std;
// Function to load graph from .txt file
Graph loadGraph(const std::string &filename, int &numNodes, int &numEdges) {
    Graph G;
    std::ifstream infile(filename);
    int u, v;

    // Read each edge and add both directions (for undirected graph)
    while (infile >> u >> v) {
        if (G[u].insert(v).second) {  // Add edge u -> v if it's new
            numEdges++;  // Count each undirected edge once
        }
        G[v].insert(u);  // Add edge v -> u (undirected)
    }

    // Set the number of nodes
    numNodes = G.size();

    return G;
}

// Función para seleccionar el pivote con el mayor número de vecinos en `cand ∪ fini`
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

// Algoritmo TTT secuencial
void TTT(const Graph &G, Set K, Set cand, Set fini, vector<Set> &maximalCliques) {
    if (cand.empty() && fini.empty()) {
        maximalCliques.push_back(K);  // Agregar K a la lista de cliques maximales
        return;
    }

    // Selección del pivote
    int pivot = Pivot(G, cand, fini);

    // Creación del conjunto `ext` como `cand - ΓG(pivot)`
    Set ext;
    for (int v : cand) {
        if (G.at(pivot).find(v) == G.at(pivot).end()) {
            ext.insert(v);
        }
    }

    // Procesar cada vértice `q` en `ext`
    for (int q : ext) {
        Set Kq = K;
        Kq.insert(q);

        // Construir `candq` y `finiq`
        Set candq, finiq;
        for (int neighbor : G.at(q)) {
            if (cand.find(neighbor) != cand.end()) {
                candq.insert(neighbor);
            }
            if (fini.find(neighbor) != fini.end()) {
                finiq.insert(neighbor);
            }
        }

        // Actualizar `cand` y `fini` para eliminar `q`
        cand.erase(q);
        fini.insert(q);

        // Llamada recursiva
        TTT(G, Kq, candq, finiq, maximalCliques);
    }
}


// Función principal
int main() {
    int numNodes = 0, numEdges = 0;
    string dataset = "CA-GrQc.txt";
    Graph G = loadGraph(dataset, numNodes, numEdges);
    //Graph G = loadGraph("com-dblp.ungraph.txt", numNodes, numEdges);
    //Graph G = loadGraph("dataset_test.txt", numNodes, numEdges);
    std::cout << "Dataset: " << dataset << std::endl;

    // Salida del número de nodos y aristas
    std::cout << "Number of nodes: " << numNodes << std::endl;
    std::cout << "Number of edges: " << numEdges << std::endl;

    Set K, cand, fini;
    for (const auto &[node, _] : G) {
        cand.insert(node);
    }

    vector<Set> maximalCliques;

    // Ejecutar TTT secuencial
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