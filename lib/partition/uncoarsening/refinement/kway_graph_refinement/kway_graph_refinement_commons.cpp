/******************************************************************************
 * kway_graph_refinement_commons.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include <omp.h>

#include "kway_graph_refinement_commons.h"

kway_graph_refinement_commons::kway_graph_refinement_commons( PartitionConfig & config ) {
        init(config);
}

kway_graph_refinement_commons::~kway_graph_refinement_commons() {
}


unsigned kway_graph_refinement_commons::getUnderlyingK() {
    return m_local_degrees.size();
}

void kway_graph_refinement_commons::init(PartitionConfig & config) {
    m_local_degrees.resize(config.k);
    for( PartitionID i = 0; i < config.k; i++) {
        m_local_degrees[i].round        = 0;
        m_local_degrees[i].local_degree = 0;
    }

    m_round = 0;//needed for the computation of internal and external degrees
}

bool kway_graph_refinement_commons::incident_to_more_than_two_partitions(graph_access & G, NodeID & node) {
    bool ret_value = false;
    PartitionID own_partition = G.getPartitionIndex(node);
    PartitionID second_partition = INVALID_PARTITION;

    forall_out_edges(G, e, node) {
                NodeID target = G.getEdgeTarget(e);
                PartitionID target_partition = G.getPartitionIndex(target);
                if(target_partition != own_partition) {
                    if(second_partition == INVALID_PARTITION) {
                        second_partition = target_partition;
                    } else if(target_partition != second_partition) {
                        ret_value = true;
                        break;
                    }
                }

            } endfor

    return ret_value;
}

bool kway_graph_refinement_commons::int_ext_degree( graph_access & G,
                                                           const NodeID & node,
                                                           const PartitionID lhs,
                                                           const PartitionID rhs,
                                                           EdgeWeight & int_degree,
                                                           EdgeWeight & ext_degree) {


    ASSERT_EQ(lhs, G.getPartitionIndex(node));

    int_degree               = 0;
    ext_degree               = 0;
    bool update_is_difficult = false;

    forall_out_edges(G, e, node) {
                NodeID target                 = G.getEdgeTarget(e);
                PartitionID targets_partition = G.getPartitionIndex(target);

                if(targets_partition == lhs) {
                    int_degree += G.getEdgeWeight(e);
                } else if(targets_partition == rhs) {
                    ext_degree += G.getEdgeWeight(e);
                }

                if(targets_partition != lhs && targets_partition != rhs) {
                    update_is_difficult = true;
                }
            } endfor

    return update_is_difficult;
}

Gain kway_graph_refinement_commons::compute_gain(graph_access &G, NodeID &node, PartitionID &max_gainer,
                                                        EdgeWeight &ext_degree) {
    //for all incident partitions compute gain
    //return max gain and max_gainer partition
    PartitionID source_partition = G.getPartitionIndex(node);
    EdgeWeight max_degree        = 0;
    max_gainer                   = INVALID_PARTITION;


    m_round++;//can become zero again ?? really?
    forall_out_edges(G, e, node) {
                NodeID target                = G.getEdgeTarget(e);
                PartitionID target_partition = G.getPartitionIndex(target);

                if(m_local_degrees[target_partition].round == m_round) {
                    m_local_degrees[target_partition].local_degree += G.getEdgeWeight(e);
                } else {
                    m_local_degrees[target_partition].local_degree = G.getEdgeWeight(e);
                    m_local_degrees[target_partition].round = m_round;
                }


                if(m_local_degrees[target_partition].local_degree >= max_degree && target_partition != source_partition) {
                    if(m_local_degrees[target_partition].local_degree > max_degree) {
                        max_degree = m_local_degrees[target_partition].local_degree;
                        max_gainer = target_partition;
                    } else {
                        //break ties randomly
                        bool accept = random_functions::nextBool();
                        if(accept) {
                            max_degree = m_local_degrees[target_partition].local_degree;
                            max_gainer = target_partition;
                        }
                    }
                }
            } endfor

    if(max_gainer != INVALID_PARTITION) {
        ext_degree = max_degree;
    } else {
        ext_degree = 0;
    }

    if(m_local_degrees[source_partition].round != m_round) {
        m_local_degrees[source_partition].local_degree = 0;
    }

    return max_degree-m_local_degrees[source_partition].local_degree;
}


Gain kway_graph_refinement_commons::gain_internal(partitionAdjMat base, partitionAdjMat now) {

    Gain beta_sb = 1.0; // discount factor
//    int tmp = beta_sb * 1;
//    std::cout<<tmp<<std::endl;
//
//    std::cout << "fuck" << beta_sb << std::endl;
//    printf("%f", beta_sb);
    while (base.getMax() == now.getMax() && base.getMax()!= 0){
        std::pair<int, int> max_ind_base = base.getMaxInd();
        std::pair<int, int> max_ind_now = now.getMaxInd();
        base[max_ind_base.first][max_ind_base.second] = 0;
        now[max_ind_now.first][max_ind_now.second] = 0;
        beta_sb *= 0.8;
    }
    return beta_sb*(base.getMax() - now.getMax());
}

Gain
kway_graph_refinement_commons::compute_gain_maxcut(graph_access &G, NodeID &node, PartitionID &max_gainer,
                                                   EdgeWeight &ext_degree, partitionAdjMat mat,
                                                   complete_boundary &bnd) {
    PartitionID source_partition = G.getPartitionIndex(node);
    Gain max_gain = -1000;
    max_gainer                   = INVALID_PARTITION;


    partitionAdjMat differceMat(G.get_partition_count());

//    m_round++;
    forall_out_edges(G, e, node) {
                NodeID target                = G.getEdgeTarget(e);
                PartitionID target_partition = G.getPartitionIndex(target);
                differceMat[source_partition][target_partition] += G.getEdgeWeight(e);
            } endfor
//    differceMat.showMat();
    // set ext_dregee
    ext_degree = 0;
    forall_blocks(G,p) {
                if (p == source_partition) continue;
                ext_degree = std::max(ext_degree, differceMat[source_partition][p]);
            } endfor

    forall_blocks(G, target){ // if we move to this partition
                if (target == source_partition) continue;

                partitionAdjMat thismat(G.get_partition_count());
                thismat.setMat(mat.getMat());
                // update two entries in the mat
                thismat[source_partition][target] += differceMat[source_partition][source_partition] - differceMat[source_partition][target];
                thismat[target][source_partition] += differceMat[source_partition][source_partition] - differceMat[source_partition][target];
                // update 2 * (n-2) entries
                forall_blocks(G,p) {
                            if (p == source_partition || p == target) continue;
                            thismat[source_partition][p] -= differceMat[source_partition][p];
                            thismat[p][source_partition] -= differceMat[source_partition][p];
                            thismat[target][p] += differceMat[source_partition][p];
                            thismat[p][target] += differceMat[source_partition][p];
                        } endfor
//        thismat.showMat();
                Gain gain_now = gain_internal(mat, thismat);
                if (max_gain < gain_now){
                    max_gain = gain_now;
                    max_gainer = target;
                }
                if (max_gain == gain_now){
                    max_gainer = bnd.getBlockWeight(target) < bnd.getBlockWeight(max_gainer) ? target : max_gainer;
                }
            } endfor
    return max_gain;
}
