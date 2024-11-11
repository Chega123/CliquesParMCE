#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <tbb/tbb.h>
#include <tbb/concurrent_vector.h>
#include <tbb/concurrent_unordered_set.h>
#include <tbb/parallel_for_each.h>
#include <dirent.h>


// Type definitions for ease
using Graph = std::unordered_map<int, std::unordered_set<int>>;
using Set = std::unordered_set<int>;
using namespace std;

// Hash function for unordered_set of integers
struct SetHash {
    size_t operator()(const Set &s) const {
        size_t hash = 0;
        for (int v : s) {
            hash ^= std::hash<int>{}(v) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        }
        return hash;
    }
};

// Function to load graph from .txt file
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

// Función para seleccionar el pivote con el mayor número de vecinos en `cand`
int ParPivot(const Graph &G, const Set &cand, const Set &fini) {
    // Unir cand y fini
    Set candFini = cand;
    candFini.insert(fini.begin(), fini.end());

    // Vector concurrente para almacenar los pares (vértice, tamaño de intersección)
    tbb::concurrent_vector<std::pair<int, int>> intersectionSizes;

    // Calculamos en paralelo el tamaño de la intersección para cada vértice
    tbb::parallel_for_each(candFini.begin(), candFini.end(), [&](int w) {
        int tw = 0;
        // Calcular la intersección entre `cand` y los vecinos de `w`
        for (int neighbor : G.at(w)) {
            if (cand.find(neighbor) != cand.end()) {
                tw++;
            }
        }
        // Guardar el vértice y el tamaño de la intersección
        intersectionSizes.push_back({w, tw});
    });

    // Encontrar el vértice con el máximo tamaño de intersección
    auto maxElement = std::max_element(intersectionSizes.begin(), intersectionSizes.end(),
        [](const std::pair<int, int> &a, const std::pair<int, int> &b) {
            return a.second < b.second;
        });

    // Retornar el vértice con la máxima intersección
    return maxElement->first;
}

// Función optimizada ParTTT para contar solo cliques maximales sin duplicados
void ParTTT(const Graph &G, Set K, Set cand, Set fini, tbb::concurrent_unordered_set<Set, SetHash> &maximalCliques) {
    if (cand.empty() && fini.empty()) {
        std::vector<int> sortedK(K.begin(), K.end());
        std::sort(sortedK.begin(), sortedK.end());
        Set sortedSet(sortedK.begin(), sortedK.end());
        maximalCliques.insert(sortedSet);
        return;
    }

    int pivot = ParPivot(G, cand,fini);
    Set ext;

    for (int v : cand) {
        if (G.at(pivot).find(v) == G.at(pivot).end()) {
            ext.insert(v);
        }
    }

    tbb::parallel_for_each(ext.begin(), ext.end(), [&](int q) {
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
int countTriangles(const Graph &G, int v) {
    int count = 0;
    for (int u : G.at(v)) {
        for (int w : G.at(u)) {
            if (G.at(v).find(w) != G.at(v).end()) {
                count++;
            }
        }
    }
    return count / 2;
}

// Función para calcular la degeneración
std::unordered_map<int, int> calculateDegeneracies(const Graph &G) {
    std::unordered_map<int, int> degeneracies;
    for (const auto &entry : G) {
        int v = entry.first;
        degeneracies[v] = G.at(v).size();
    }
    return degeneracies;
}

// Función para computar el ranking
std::unordered_map<int, std::pair<int, int>> computeRank(const Graph &G, const std::string &method) {
    std::unordered_map<int, std::pair<int, int>> ranks;

    if (method == "degree") {
        for (const std::pair<int, std::unordered_set<int>> &entry : G) {
            int v = entry.first;
            int degree = G.at(v).size();
            ranks[v] = std::make_pair(degree, v);
        }
    } else if (method == "triangle") {
        for (const std::pair<int, std::unordered_set<int>> &entry : G) {
            int v = entry.first;
            int triangle_count = countTriangles(G, v);
            ranks[v] = std::make_pair(triangle_count, v);
        }
    } else if (method == "degeneracy") {
        std::unordered_map<int, int> degeneracies = calculateDegeneracies(G);
        for (const std::pair<int, std::unordered_set<int>> &entry : G) {
            int v = entry.first;
            ranks[v] = std::make_pair(degeneracies[v], v);
        }
    }

    return ranks;
}

void ParMCE(const Graph &G, std::unordered_map<int, std::pair<int, int>> &ranks,
            tbb::concurrent_unordered_set<Set, SetHash> &maximalCliques) {
    tbb::parallel_for_each(G.begin(), G.end(), [&](const std::pair<int, std::unordered_set<int>> &entry) {
        int v = entry.first;

        // Crear el subgrafo inducido por v y su vecindad
        Set GvNeighbors = G.at(v);
        GvNeighbors.insert(v);

        Set K = {v};
        Set cand, fini;

        for (int w : GvNeighbors) {
            if (ranks[w].first >= ranks[v].first) {
                cand.insert(w);
            } else {
                fini.insert(w);
            }
        }

        // Llamada a ParTTT para enumerar los cliques maximales
        ParTTT(G, K, cand, fini, maximalCliques);
    });
}


int main() {
    std::string dataset_folder = "/home/chega/cliquesparalela/dataset";
    std::vector<int> thread_counts = {2, 4, 8, 16};
    std::vector<std::string> methods = {"degeneracy", "degree", "triangle"};
    
    // Abrir el directorio usando dirent.h
    DIR *dir = opendir(dataset_folder.c_str());
    if (dir == nullptr) {
        std::cerr << "Error opening directory" << std::endl;
        return 1;
    }

    struct dirent *entry;
    while ((entry = readdir(dir)) != nullptr) {
        // Ignorar directorios "." y ".."
        if (entry->d_type == DT_DIR) continue;

        // Crear la ruta completa al archivo
        std::string dataset = dataset_folder + "/" + entry->d_name;
        std::cout << "Processing dataset: " << dataset << std::endl;

        int numNodes = 0, numEdges = 0;
        Graph G = loadGraph(dataset, numNodes, numEdges);
        std::cout << "Number of nodes: " << numNodes << ", Number of edges: " << numEdges << std::endl;
                                                                                
        for (const auto& method : methods) {
            std::cout << "Method: " << method << std::endl;
            auto ranks = computeRank(G, method);

            for (int num_threads : thread_counts) {
                tbb::global_control control(tbb::global_control::max_allowed_parallelism, num_threads);
                tbb::concurrent_unordered_set<Set, SetHash> maximalCliques;

                tbb::tick_count start = tbb::tick_count::now();
                ParMCE(G, ranks, maximalCliques);
                tbb::tick_count end = tbb::tick_count::now();
                double elapsedTime = (end - start).seconds();

                std::cout << "Threads: " << num_threads << ", Time: " << elapsedTime << " seconds, ";
                std::cout << "Total maximal cliques: " << maximalCliques.size() << std::endl;
            }
        }
    }

    closedir(dir);  // Cerrar el directorio después de recorrerlo
    return 0;
}