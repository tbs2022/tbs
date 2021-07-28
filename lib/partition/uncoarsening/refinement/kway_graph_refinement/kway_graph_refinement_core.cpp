/******************************************************************************
 * kway_graph_refinement_core.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include <algorithm>

#include "data_structure/priority_queues/bucket_pq.h"
#include "data_structure/priority_queues/maxNodeHeap.h"
#include "kway_graph_refinement_core.h"
#include "kway_stop_rule.h"
#include "quality_metrics.h"
#include "random_functions.h"

kway_graph_refinement_core::kway_graph_refinement_core() : commons (NULL){
}

kway_graph_refinement_core::~kway_graph_refinement_core() {
        if( commons != NULL)  delete commons;
}
EdgeWeight kway_graph_refinement_core::single_kway_refinement_round(PartitionConfig & config, 
                                                                    graph_access & G, 
                                                                    complete_boundary & boundary, 
                                                                    boundary_starting_nodes & start_nodes, 
                                                                    int step_limit, 
                                                                    vertex_moved_hashtable & moved_idx) {
        std::unordered_map<PartitionID, PartitionID> touched_blocks;
        return single_kway_refinement_round_internal(config, G, boundary, start_nodes, 
                                                     step_limit, moved_idx, false, touched_blocks);
}

EdgeWeight kway_graph_refinement_core::single_kway_refinement_round(PartitionConfig & config, 
                                                                    graph_access & G, 
                                                                    complete_boundary & boundary, 
                                                                    boundary_starting_nodes & start_nodes, 
                                                                    int step_limit, 
                                                                    vertex_moved_hashtable & moved_idx,
                                                                    std::unordered_map<PartitionID, PartitionID> & touched_blocks) {

        return single_kway_refinement_round_internal(config, G, boundary, start_nodes, 
                                                     step_limit, moved_idx, true, touched_blocks);
}


EdgeWeight kway_graph_refinement_core::single_kway_refinement_round_internal(PartitionConfig & config, 
                                                                    graph_access & G, 
                                                                    complete_boundary & boundary, 
                                                                    boundary_starting_nodes & start_nodes, 
                                                                    int step_limit,
                                                                    vertex_moved_hashtable & moved_idx,
                                                                    bool compute_touched_partitions,
                                                                    std::unordered_map<PartitionID, PartitionID> &  touched_blocks) {

        if( commons == NULL ) commons = new kway_graph_refinement_commons(config);

        refinement_pq* queue = NULL;
        if(config.use_bucket_queues) {
                EdgeWeight max_degree = G.getMaxDegree();
                queue                 = new bucket_pq(max_degree);
        } else {
                queue                 = new maxNodeHeap(); 
        }

        init_queue_with_boundary(config, G, start_nodes, queue, moved_idx, boundary);
        std::cout << "Queue size = " << queue->size() << std::endl;
        
        if(queue->empty()) {delete queue; return 0;}

        std::vector<NodeID> transpositions;
        std::vector<PartitionID> from_partitions;
        std::vector<PartitionID> to_partitions;

        int max_number_of_swaps = (int)(G.number_of_nodes());
        int min_cut_index       = -1;
        quality_metrics qm;
        Gain cut         = qm.max_pairwise_cut(G); // so we dont need to compute the edge cut
        Gain initial_cut = cut;

        //roll forwards
        Gain best_cut = cut;
        int number_of_swaps = 0;
        int movements       = 0;

        kway_stop_rule* stopping_rule = NULL;
        switch(config.kway_stop_rule) {
                case KWAY_SIMPLE_STOP_RULE: 
                        stopping_rule = new kway_simple_stop_rule(config);
                        break;
//                case KWAY_ADAPTIVE_STOP_RULE:
//                        stopping_rule = new kway_adaptive_stop_rule(config);
//                        break;

        }

        for(number_of_swaps = 0, movements = 0; movements < max_number_of_swaps; movements++, number_of_swaps++) {
                if( queue->empty() ) {std::cout << "Queue empty" << std::endl; break;}
//                std::cout << number_of_swaps - min_cut_index <<std::endl;
                if( stopping_rule->search_should_stop(min_cut_index, number_of_swaps, step_limit) ) break;
                Gain gain = queue->maxValue();
                NodeID node = queue->deleteMax();
//                if (gain <= 0) { std::cout << "No positive gain" << std::endl; break;} // change here to control what kind of nodes can move

#ifndef NDEBUG
                //a new partition adj matrix every iteration
                partitionAdjMat mat(G.get_partition_count());
                mat.pairwise_cut(G);
//                mat.showMat();
                PartitionID maxgainer;
                EdgeWeight ext_degree;
                ASSERT_TRUE(moved_idx[node].index == NOT_MOVED);
                ASSERT_EQ(gain, commons->compute_gain_maxcut(G, node, maxgainer, ext_degree,mat, boundary));
                ASSERT_TRUE(ext_degree > 0);
//                assert(gain_statisfy(gain));
#endif

                PartitionID from = G.getPartitionIndex(node); 
                bool successfull = move_node(config, G, node, moved_idx, queue, boundary);

                if(successfull) {
                        cut = qm.max_pairwise_cut(G);
                        stopping_rule->push_statistics(gain);

                        bool accept_equal = random_functions::nextBool();
accept_equal = false;
                        if( cut < best_cut || ( cut == best_cut && accept_equal )) {
//                               std::cout << "cut now = " << cut << std::endl;
                                best_cut = cut;
                                min_cut_index = number_of_swaps;
                                if(cut < best_cut)
                                        stopping_rule->reset_statistics();
                        }

                        from_partitions.push_back(from);
                        to_partitions.push_back(G.getPartitionIndex(node));
                        transpositions.push_back(node);
                } else {
                        number_of_swaps--; //because it wasnt swaps
                }
                moved_idx[node].index = MOVED;
                ASSERT_TRUE(boundary.assert_bnodes_in_boundaries());
                ASSERT_TRUE(boundary.assert_boundaries_are_bnodes());

        } 

        ASSERT_TRUE(boundary.assert_bnodes_in_boundaries());
        ASSERT_TRUE(boundary.assert_boundaries_are_bnodes());
/////////////////////////////////////////
/*   no need to roll back if gain > 0
        //roll backwards
        for(number_of_swaps--; number_of_swaps>min_cut_index; number_of_swaps--) {
                ASSERT_TRUE(transpositions.size() > 0);

                NodeID node = transpositions.back();
                transpositions.pop_back();

                PartitionID to = from_partitions.back();
                from_partitions.pop_back();
                to_partitions.pop_back();

                move_node_back(config, G, node, to,  moved_idx, queue, boundary);
        }
*/
       
        //reconstruct the touched partitions
        if(compute_touched_partitions) {
                ASSERT_EQ(from_partitions.size(), to_partitions.size());
                for(unsigned i = 0; i < from_partitions.size(); i++) {
                        touched_blocks[from_partitions[i]] = from_partitions[i];
                        touched_blocks[to_partitions[i]]   = to_partitions[i];
                }
        }

        ASSERT_TRUE(boundary.assert_bnodes_in_boundaries());
        ASSERT_TRUE(boundary.assert_boundaries_are_bnodes());

    //a new partition adj matrix every iteration
    partitionAdjMat mat(G.get_partition_count());
    mat.pairwise_cut(G);
    std::cout << "max cut" << mat.getMax() << std::endl;
    delete queue;
    delete stopping_rule;
        return initial_cut - best_cut; 
}

void kway_graph_refinement_core::init_queue_with_boundary(const PartitionConfig &config, graph_access &G,
                                                          std::vector<NodeID> &bnd_nodes, refinement_pq *queue,
                                                          vertex_moved_hashtable &moved_idx, complete_boundary &bnd) {

        if(config.permutation_during_refinement == PERMUTATION_QUALITY_FAST) {
                random_functions::permutate_vector_fast(bnd_nodes, false);
        } else if(config.permutation_during_refinement == PERMUTATION_QUALITY_GOOD) {
                random_functions::permutate_vector_good(bnd_nodes, false);
        }

        partitionAdjMat mat(G.get_partition_count());
        mat.pairwise_cut(G);
//        mat.showMat();

        for( unsigned int i = 0; i < bnd_nodes.size(); i++) {
                NodeID node = bnd_nodes[i];

                if( moved_idx.find(node) == moved_idx.end() ) {
                        PartitionID max_gainer;
                        EdgeWeight ext_degree;
                        //compute gain
//                        Gain gain = commons->compute_gain(G, node, max_gainer, ext_degree);
                    Gain gain = commons->compute_gain_maxcut(G, node, max_gainer, ext_degree, mat, bnd);
                    if (gain_statisfy(gain)) {
                        queue->insert(node, gain);
                        moved_idx[node].index = NOT_MOVED;
                    }
                }
        }
}


void kway_graph_refinement_core::move_node_back(PartitionConfig & config, 
                graph_access & G, 
                NodeID & node,
                PartitionID & to, 
                vertex_moved_hashtable & moved_idx, 
                refinement_pq * queue, 
                complete_boundary & boundary) {

        PartitionID from = G.getPartitionIndex(node);
        G.setPartitionIndex(node, to);        

        boundary_pair pair;
        pair.k   = config.k;
        pair.lhs = from;
        pair.rhs = to;

        //update all boundaries
        boundary.postMovedBoundaryNodeUpdates(node, &pair, true, true);

        NodeWeight this_nodes_weight = G.getNodeWeight(node);
        boundary.setBlockNoNodes(from, boundary.getBlockNoNodes(from)-1);
        boundary.setBlockNoNodes(to,   boundary.getBlockNoNodes(to)+1);
        boundary.setBlockWeight( from, boundary.getBlockWeight(from)-this_nodes_weight);
        boundary.setBlockWeight( to,   boundary.getBlockWeight(to)+this_nodes_weight);
}

void kway_graph_refinement_core::setup_start_nodes(PartitionConfig & config, graph_access & G, complete_boundary & boundary, boundary_starting_nodes & start_nodes) {
    QuotientGraphEdges quotient_graph_edges;
    boundary.getQuotientGraphEdges(quotient_graph_edges);

    std::unordered_map<NodeID, bool> allready_contained;

    for( unsigned i = 0; i < quotient_graph_edges.size(); i++) {
        boundary_pair & ret_value = quotient_graph_edges[i];
        PartitionID lhs           = ret_value.lhs;
        PartitionID rhs           = ret_value.rhs;

        PartialBoundary & partial_boundary_lhs = boundary.getDirectedBoundary(lhs, lhs, rhs);
        forall_boundary_nodes(partial_boundary_lhs, cur_bnd_node) {
                    ASSERT_EQ(G.getPartitionIndex(cur_bnd_node), lhs);
                    if(allready_contained.find(cur_bnd_node) == allready_contained.end() ) {
                        start_nodes.push_back(cur_bnd_node);
                        allready_contained[cur_bnd_node] = true;
                    }
                } endfor

        PartialBoundary & partial_boundary_rhs = boundary.getDirectedBoundary(rhs, lhs, rhs);
        forall_boundary_nodes(partial_boundary_rhs, cur_bnd_node) {
                    ASSERT_EQ(G.getPartitionIndex(cur_bnd_node), rhs);
                    if(allready_contained.find(cur_bnd_node) == allready_contained.end()) {
                        start_nodes.push_back(cur_bnd_node);
                        allready_contained[cur_bnd_node] = true;
                    }
                } endfor
    }
}
inline bool kway_graph_refinement_core::move_node(PartitionConfig & config,
                                                  graph_access & G,
                                                  NodeID & node,
                                                  vertex_moved_hashtable & moved_idx,
                                                  refinement_pq * queue,
                                                  complete_boundary & boundary) {

    ASSERT_TRUE(boundary.assert_bnodes_in_boundaries());
    ASSERT_TRUE(boundary.assert_boundaries_are_bnodes());

    PartitionID from = G.getPartitionIndex(node);
    PartitionID to;
    EdgeWeight node_ext_deg;
    partitionAdjMat mat0(G.get_partition_count());
    mat0.pairwise_cut(G);
    commons->compute_gain_maxcut(G, node, to, node_ext_deg,mat0,boundary); //fill  partitioin ID to

    NodeWeight this_nodes_weight = G.getNodeWeight(node);
//    std::cout << "upper bound is in move nodes>" << config.upper_bound_partition << std::endl;
//    if(boundary.getBlockWeight(to) + this_nodes_weight >= config.upper_bound_partition)
//        return false;

    if(boundary.getBlockWeight(to) + this_nodes_weight >= 1000)
        return false;

    if(boundary.getBlockNoNodes(from) - 1 == 0) // assure that no block gets accidentally empty
        return false;

    /* G changes here
     *
     *
     *
     *
     * */
    G.setPartitionIndex(node, to);

    boundary_pair pair;
    pair.k = config.k;
    pair.lhs = from;
    pair.rhs = to;

    /* boundary changes here
     *
     *
     *
     *
     * */
    boundary.postMovedBoundaryNodeUpdates(node, &pair, true, true);

    boundary.setBlockNoNodes(from, boundary.getBlockNoNodes(from)-1);
    boundary.setBlockNoNodes(to,   boundary.getBlockNoNodes(to)+1);
    boundary.setBlockWeight( from, boundary.getBlockWeight(from)-this_nodes_weight);
    boundary.setBlockWeight( to,   boundary.getBlockWeight(to)+this_nodes_weight);

    ASSERT_TRUE(boundary.assert_bnodes_in_boundaries());
    ASSERT_TRUE(boundary.assert_boundaries_are_bnodes());

    // G is already changed so a new partition adj matrix is necessary
    partitionAdjMat mat(G.get_partition_count());
    mat.pairwise_cut(G);

//    std::cout<< mat0.getMax() << " " << mat.getMax() << std::endl;
//    ASSERT_GEQ(mat0.getMax(), mat.getMax());
    // set up all boundary nodes
    boundary_starting_nodes all_boundary_nodes;
    setup_start_nodes(config,G,boundary,all_boundary_nodes);
    // all boundary needs update
    for (auto target : all_boundary_nodes) {
        PartitionID targets_max_gainer;
        EdgeWeight ext_degree; // the local external degree
        Gain gain = commons->compute_gain_maxcut(G, target, targets_max_gainer, ext_degree, mat, boundary);

        if(queue->contains(target)) {
            assert(moved_idx.find(target) != moved_idx.end());
            if(ext_degree > 0) {
                queue->changeKey(target, gain);
            } else {
                queue->deleteNode(target);
            }
        } else {
            /////////////////////////////////////////////////////////attaintion
            if(ext_degree > 0 && gain_statisfy(gain)) {
                if(moved_idx.find(target) == moved_idx.end()) { // if not moved before
                    queue->insert(target, gain);
                    moved_idx[target].index = NOT_MOVED;
                }
            }
        }
    }
    //update gain of neighbors / the boundaries have allready been updated

    forall_out_edges(G, e, node) {
                NodeID target = G.getEdgeTarget(e);
                PartitionID targets_max_gainer;
                EdgeWeight ext_degree; // the local external degree
//                Gain gain = commons->compute_gain(G, target, targets_max_gainer, ext_degree);
                Gain gain = commons->compute_gain_maxcut(G, target, targets_max_gainer, ext_degree, mat, boundary);

                if(queue->contains(target)) {
                    assert(moved_idx.find(target) != moved_idx.end());
                    if(ext_degree > 0) {
                        queue->changeKey(target, gain);
                    } else {
                        queue->deleteNode(target);
                    }
                } else {
                    /////////////////////////////////////////////////////////attaintion
                    if(ext_degree > 0 && gain_statisfy(gain)){
                        if(moved_idx.find(target) == moved_idx.end()) { // if not moved before
                            queue->insert(target, gain);
                            moved_idx[target].index = NOT_MOVED;
                        }
                    }
                }
            } endfor

    return true;
}

bool kway_graph_refinement_core::gain_statisfy(Gain gain) {
    return gain >= 0;
}
