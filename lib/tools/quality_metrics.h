/******************************************************************************
 * quality_metrics.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/


#ifndef QUALITY_METRICS_10HC2I5M
#define QUALITY_METRICS_10HC2I5M

#include "data_structure/graph_access.h"
#include "data_structure/matrix/matrix.h"
#include "partition_config.h"
#define Matrix std::vector<std::vector<EdgeWeight>>
class partitionAdjMat {
public:
    int n; // adjMat is n x n

    partitionAdjMat(int n) : n(n){
        adjMat.resize(n);
        for (int i = 0; i < n; ++i) {
            adjMat[i] = std::vector<EdgeWeight>(n, 0);
        };
    }

    void showMat() {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                std::cout << adjMat[i][j] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    Matrix getMat() {
        return adjMat;
    }
    void setMat(Matrix mat) {
        adjMat = mat;
    }
    std::vector<int> &operator[] (int i) {
        adjMat[i];
    }
    EdgeWeight getMax() {
        EdgeWeight maxnum = 0;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                maxnum = std::max(maxnum, adjMat[i][j]);
            }
        }
        return maxnum;
    }
    std::pair<int,int> getMaxInd() {
        EdgeWeight maxnum = 0;
        std::pair<int, int> index;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (maxnum < adjMat[i][j]){
                    index.first = i;
                    index.second = j;
                }
            }
        }
        return index;
    }
    void pairwise_cut(graph_access &G);

private:
    Matrix adjMat;
};

class quality_metrics {
public:
        quality_metrics();
        virtual ~quality_metrics ();

        EdgeWeight edge_cut(graph_access & G);
        EdgeWeight edge_cut(graph_access & G, int * partition_map); 
        EdgeWeight edge_cut(graph_access & G, PartitionID lhs, PartitionID rhs);
        EdgeWeight max_communication_volume(graph_access & G);
        EdgeWeight min_communication_volume(graph_access & G);
        EdgeWeight max_communication_volume(graph_access & G, int * partition_map);
        EdgeWeight total_communication_volume(graph_access & G); 
        EdgeWeight objective(const PartitionConfig & config, graph_access & G, int * partition_map);
        EdgeWeight edge_cut_connected(graph_access & G, int * partition_map);
        int boundary_nodes(graph_access & G);
        NodeWeight separator_weight(graph_access& G);
        double balance(graph_access & G);
        double balance_edges(graph_access & G);
        double balance_separator(graph_access & G);
        double edge_balance(graph_access &G, const std::vector<PartitionID> &edge_partition);

        NodeWeight total_qap(graph_access & C, matrix & D, std::vector< NodeID > & rank_assign);
        NodeWeight total_qap(matrix & C, matrix & D, std::vector< NodeID > & rank_assign);

    EdgeWeight max_pairwise_cut(graph_access &G);

    void pairwise_cut(graph_access &G, partitionAdjMat thisMat);
};

#endif /* end of include guard: QUALITY_METRICS_10HC2I5M */
