#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <tbb/tbb.h>
#include <tbb/concurrent_vector.h>
#include <limits>
#include <tbb/blocked_range.h>
#include <tbb/concurrent_unordered_set.h>
#include <tbb/parallel_for_each.h>
#include <dirent.h>
#include <atomic>
#include <algorithm>

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
    if (G.find(v) == G.end()) {
        return 0;
    }
    int count = 0;
    const auto &neighbors = G.at(v);
    std::vector<int> local_counts(neighbors.size(), 0);
    std::vector<int> neighbor_vector(neighbors.begin(), neighbors.end());
    tbb::parallel_for(tbb::blocked_range<size_t>(0, neighbor_vector.size()), [&](const tbb::blocked_range<size_t> &range) {
        for (size_t i = range.begin(); i < range.end(); ++i) {
            int u = neighbor_vector[i];
            int local_count = 0;

            if (G.find(u) != G.end()) {
                for (int w : G.at(u)) {
                    if (neighbors.find(w) != neighbors.end()) {
                        local_count++;
                    }
                }
            }

            local_counts[i] = local_count; 
        }
    });
    for (int local_count : local_counts) {
        count += local_count;
    }

    return count / 2; 
}

int calculateDegeneracyTBB(const Graph &G, int node) {
    if (G.find(node) == G.end()) {
        return 0; 
    }
    int min_degree = G.at(node).size(); 
    std::vector<int> neighbors(G.at(node).begin(), G.at(node).end());
    std::atomic<int> local_min(min_degree); 

    tbb::parallel_for(tbb::blocked_range<size_t>(0, neighbors.size()), 
        [&](const tbb::blocked_range<size_t>& range) {
            int local_min_value = local_min.load();  
            for (size_t i = range.begin(); i < range.end(); ++i) {
                int neighbor = neighbors[i];
                if (G.find(neighbor) != G.end()) {
                    int neighbor_degree = G.at(neighbor).size();
                    local_min_value = std::min(local_min_value, neighbor_degree);
                }
            }
            local_min.compare_exchange_strong(local_min_value, std::min(local_min_value, local_min.load()));
        }
    );

    return local_min.load();
}

std::unordered_map<int, std::pair<int, int>> computeRank(const Graph &G, const std::string &method) {
    tbb::concurrent_unordered_map<int, std::pair<int, int>> concurrent_ranks;

    if (method == "degree") {
        tbb::parallel_for_each(G.begin(), G.end(), [&](const std::pair<const int, std::unordered_set<int>> &entry) {
            int v = entry.first;
            int degree = G.at(v).size();
            concurrent_ranks[v] = std::make_pair(degree, v);
        });

    } else if (method == "triangle") {
        tbb::parallel_for_each(G.begin(), G.end(), [&](const std::pair<const int, std::unordered_set<int>> &entry) {
            int v = entry.first;
            int triangle_count = countTriangles(G, v);
            concurrent_ranks[v] = std::make_pair(triangle_count, v);
        });
    }else if (method == "degeneracy") {
        tbb::parallel_for_each(G.begin(), G.end(), [&](const std::pair<const int, std::unordered_set<int>> &entry) {
            int v = entry.first;
            int degeneracy = calculateDegeneracyTBB(G, v); 
            concurrent_ranks[v] = std::make_pair(degeneracy, v);
        });
    }
    // Convertir concurrent_unordered_map a unordered_map
    std::unordered_map<int, std::pair<int, int>> ranks(concurrent_ranks.begin(), concurrent_ranks.end());
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
    std::string dataset_folder = "/home/chega/cliquesparalela/datasets_prubas";
    std::vector<int> thread_counts = {2,4,8,16};
    std::vector<std::string> methods = { "degree", "triangle","degeneracy"};
    
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