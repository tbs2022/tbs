/******************************************************************************
 * kway_graph_refinement_commons.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef KWAY_GRAPH_REFINEMENT_COMMONS_PVGY97EW
#define KWAY_GRAPH_REFINEMENT_COMMONS_PVGY97EW

#include <vector>

#include "data_structure/priority_queues/priority_queue_interface.h"
#include "definitions.h"
#include "random_functions.h"
#include "uncoarsening/refinement/refinement.h"
#include "uncoarsening/refinement/quotient_graph_refinement/2way_fm_refinement/vertex_moved_hashtable.h"
#include "quality_metrics.h"

class kway_graph_refinement_commons  {
        public:

                kway_graph_refinement_commons( PartitionConfig & config );
                virtual ~kway_graph_refinement_commons();

                void init( PartitionConfig & config );

                bool incident_to_more_than_two_partitions(graph_access & G, NodeID & node);

                Gain compute_gain(graph_access &G, NodeID &node, PartitionID &max_gainer, EdgeWeight &ext_degree);

                Gain
                compute_gain_maxcut(graph_access &G, NodeID &node, PartitionID &max_gainer,
                                    EdgeWeight &ext_degree, partitionAdjMat mat,
                                    complete_boundary &bnd);

                bool int_ext_degree( graph_access & G, 
                                     const NodeID & node,
                                     const PartitionID lhs,
                                     const PartitionID rhs,
                                     EdgeWeight & int_degree,
                                     EdgeWeight & ext_degree);

                inline unsigned getUnderlyingK();

        private:

                //for efficient computation of internal and external degrees
                struct round_struct {
                        unsigned round;
                        EdgeWeight local_degree;
                };

                std::vector<round_struct>                    m_local_degrees;
                unsigned                                     m_round;

    Gain gain_internal(partitionAdjMat base, partitionAdjMat now);
};

#endif /* end of include guard: KWAY_GRAPH_REFINEMENT_COMMONS_PVGY97EW */

