#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <tbb/tbb.h>
#include <tbb/concurrent_vector.h>
#include <tbb/blocked_range.h>
#include <tbb/concurrent_unordered_set.h>
#include <tbb/parallel_for_each.h>
#include <dirent.h>

using Graph = std::unordered_map<int, std::unordered_set<int>>;
using Set = std::unordered_set<int>;
using namespace std;

// hash function para los ints en set
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

// para seleccionar el pivote con el mayor número de vecinos en cand
int ParPivot(const Graph &G, const Set &cand, const Set &fini) {
    // unir cand y fini
    Set candFini = cand;
    candFini.insert(fini.begin(), fini.end());

    // vector concurrente para almacenar los pares (vertice, tamaño de interseccion)
    tbb::concurrent_vector<std::pair<int, int>> intersectionSizes;

    // calculamos en paralelo el tamaño de la interseccion para cada vertice
    tbb::parallel_for_each(candFini.begin(), candFini.end(), [&](int w) {
        int tw = 0;
        // calcular la intersección entre cand y los vecinos de w
        for (int neighbor : G.at(w)) {
            if (cand.find(neighbor) != cand.end()) {
                tw++;
            }
        }
        // guardar el vertice y el tamaño de la interseccion
        intersectionSizes.push_back({w, tw});
    });

    // busca el vertice con el maximo tamaño de interseccion
    auto maxElement = std::max_element(intersectionSizes.begin(), intersectionSizes.end(),
        [](const std::pair<int, int> &a, const std::pair<int, int> &b) {
            return a.second < b.second;
        });

    // retorna el vertice con el maximo tamaño de interseccion
    return maxElement->first;
}


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
    
    // Usamos un vector para almacenar los conteos locales por hilo
    std::vector<int> local_counts(G.at(v).size(), 0);

    // Paralelización del conteo de triángulos
    tbb::parallel_for(tbb::blocked_range<size_t>(0, G.at(v).size()), [&](const tbb::blocked_range<size_t>& range) {
        for (size_t i = range.begin(); i < range.end(); ++i) {
            int u = *std::next(G.at(v).begin(), i); // Obtener el nodo correspondiente
            int local_count = 0;
            for (int w : G.at(u)) {
                if (G.at(v).find(w) != G.at(v).end()) {
                    local_count++;
                }
            }
            local_counts[i] = local_count; // Almacenar el conteo local
        }
    });

    // Ahora sumamos los resultados locales para obtener el conteo total
    for (int local_count : local_counts) {
        count += local_count;
    }

    return count / 2; // Dividimos entre 2 para no contar los triángulos dos veces
}




// calcular la degeneración
std::unordered_map<int, int> calculateDegeneracies(const Graph &G) {
    std::unordered_map<int, int> degeneracies;

    // Usar tbb::concurrent_unordered_map para ser seguro en paralelo
    tbb::concurrent_unordered_map<int, int> concurrent_degeneracies;

    tbb::parallel_for_each(G.begin(), G.end(), [&](const std::pair<int, std::unordered_set<int>> &entry) {
        int v = entry.first;
        int degree = G.at(v).size();
        concurrent_degeneracies[v] = degree;
    });

    // Copiar resultados a un std::unordered_map si se necesita
    for (const auto &entry : concurrent_degeneracies) {
        degeneracies[entry.first] = entry.second;
    }
    return degeneracies;
}


// ranking
std::unordered_map<int, std::pair<int, int>> computeRank(const Graph &G, const std::string &method) {
    std::unordered_map<int, std::pair<int, int>> ranks;

    if (method == "degree") {
        for (const std::pair<int, std::unordered_set<int>> &entry : G) {
            int v = entry.first;
            int degree = G.at(v).size();
            ranks[v] = std::make_pair(degree, v);
        }
    }
    else if (method == "triangle") {
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
    else if (method == "degree_degeneracy") {
        double degree_weight = 0.7;  // Peso para el grado
        double degeneracy_weight = 0.3;  // Peso para el índice de degeneración

        std::unordered_map<int, int> degeneracies = calculateDegeneracies(G);  // Calculamos la degeneración de cada nodo

        for (const std::pair<int, std::unordered_set<int>> &entry : G) {
            int v = entry.first;
            
            // Grado del nodo
            int degree = G.at(v).size();
            
            // Índice de degeneración del nodo
            int degeneracy = degeneracies[v];
            
            // Ranking ponderado combinando los dos
            double rank = degree_weight * degree + degeneracy_weight * degeneracy;
            
            // Guardamos el ranking en el mapa
            ranks[v] = std::make_pair(rank, v);
        }

    }
    return ranks;
}

void ParMCE(const Graph &G, std::unordered_map<int, std::pair<int, int>> &ranks,
            tbb::concurrent_unordered_set<Set, SetHash> &maximalCliques) {
    tbb::parallel_for_each(G.begin(), G.end(), [&](const std::pair<int, std::unordered_set<int>> &entry) {
        int v = entry.first;

        // subgrafo inducido por v y su vecindad
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

        ParTTT(G, K, cand, fini, maximalCliques);
    });
}


int main() {
    std::string dataset_folder = "/Users/gabrielrodriguezpostigo/Library/Mobile Documents/com~apple~CloudDocs/Programación/CliquesParMCE-main-2/datasets_prubas";
    std::vector<int> thread_counts = {2, 4, 8};
    std::vector<std::string> methods = {"degeneracy", "degree", "triangle","degree_degeneracy"};
    
    DIR *dir = opendir(dataset_folder.c_str());
    if (dir == nullptr) {
        std::cerr << "Error opening directory" << std::endl;
        return 1;
    }

    struct dirent *entry;
    while ((entry = readdir(dir)) != nullptr) {
        if (entry->d_type == DT_DIR) continue;
        std::string dataset = dataset_folder + "/" + entry->d_name;
        std::cout << "Processing dataset: " << dataset << std::endl;

        int numNodes = 0, numEdges = 0;
        Graph G = loadGraph(dataset, numNodes, numEdges);
        std::cout << "Number of nodes: " << numNodes << ", Number of edges: " << numEdges << std::endl;
                                                                                
        for (const auto& method : methods) {
            std::cout << "Method: " << method << std::endl;
            

            for (int num_threads : thread_counts) {
                tbb::global_control control(tbb::global_control::max_allowed_parallelism, num_threads);
                tbb::concurrent_unordered_set<Set, SetHash> maximalCliques;

                tbb::tick_count start = tbb::tick_count::now();
                auto ranks = computeRank(G, method);
                ParMCE(G, ranks, maximalCliques);
                tbb::tick_count end = tbb::tick_count::now();
                double elapsedTime = (end - start).seconds();

                std::cout << "Threads: " << num_threads << ", Time: " << elapsedTime << " seconds, ";
                std::cout << "Total maximal cliques: " << maximalCliques.size() << std::endl;
            }
        }
    }

    closedir(dir);
    return 0;
}