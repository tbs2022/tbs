/******************************************************************************
 * kway_graph_refinement_core.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef KWAY_GRAPH_REFINEMENT_CORE_PVGY97EW
#define KWAY_GRAPH_REFINEMENT_CORE_PVGY97EW

#include <unordered_map>
#include <vector>

#include "data_structure/priority_queues/priority_queue_interface.h"
#include "definitions.h"
#include "kway_graph_refinement_commons.h"
#include "tools/random_functions.h"
#include "uncoarsening/refinement/quotient_graph_refinement/2way_fm_refinement/vertex_moved_hashtable.h"
#include "uncoarsening/refinement/refinement.h"

class kway_graph_refinement_core  {
        public:
                kway_graph_refinement_core( );
                virtual ~kway_graph_refinement_core();

                EdgeWeight single_kway_refinement_round(PartitionConfig & config, 
                                                        graph_access & G, 
                                                        complete_boundary & boundary, 
                                                        boundary_starting_nodes & start_nodes, 
                                                        int step_limit, 
                                                        vertex_moved_hashtable & moved_idx );

                EdgeWeight single_kway_refinement_round(PartitionConfig & config, 
                                                        graph_access & G, 
                                                        complete_boundary & boundary, 
                                                        boundary_starting_nodes & start_nodes, 
                                                        int step_limit, 
                                                        vertex_moved_hashtable & moved_idx,
                                                        std::unordered_map<PartitionID, PartitionID> & touched_blocks); 


         private:
               EdgeWeight single_kway_refinement_round_internal(PartitionConfig & config, 
                                                                graph_access & G, 
                                                                complete_boundary & boundary, 
                                                                boundary_starting_nodes & start_nodes, 
                                                                int step_limit,
                                                                vertex_moved_hashtable & moved_idx,
                                                                bool compute_touched_partitions,
                                                                std::unordered_map<PartitionID, PartitionID> &  touched_blocks); 


                void init_queue_with_boundary(const PartitionConfig &config, graph_access &G,
                                              std::vector<NodeID> &bnd_nodes, refinement_pq *queue,
                                              vertex_moved_hashtable &moved_idx, complete_boundary &bnd);

                inline bool move_node(PartitionConfig & config, 
                                      graph_access & G, 
                                      NodeID & node, 
                                      vertex_moved_hashtable & moved_idx, 
                                      refinement_pq * queue, 
                                      complete_boundary & boundary);

                inline void move_node_back(PartitionConfig & config, 
                                           graph_access & G, 
                                           NodeID & node,
                                           PartitionID & to, 
                                           vertex_moved_hashtable & moved_idx, 
                                           refinement_pq * queue, 
                                           complete_boundary & boundary);

                void initialize_partition_moves_array(PartitionConfig & config, 
                                                      complete_boundary & boundary, 
                                                      std::vector<bool> & partition_move_valid);

                kway_graph_refinement_commons* commons;

    void setup_start_nodes(PartitionConfig &config, graph_access &G, complete_boundary &boundary,
                           boundary_starting_nodes &start_nodes);

    bool gain_statisfy(Gain gain);
};




#endif //[> end of include guard: KWAY_GRAPH_REFINEMENT_PVGY97EW <]

