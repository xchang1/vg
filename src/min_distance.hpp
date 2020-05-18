#ifndef VG_MIN_DISTANCE_HPP_INCLUDED
#define VG_MIN_DISTANCE_HPP_INCLUDED

#include <unordered_set>

#include "snarls.hpp"
#include "hash_map.hpp"

#include "bdsg/hash_graph.hpp"

using namespace sdsl;
namespace vg { 

/**
 * The distance index. Stores minimum distances among nodes in each 
 * netgraph and chain. Used for calculation of the minimum distance between
 * two positions and for a maximum distance estimation. The maximum distance
 * estimation is at least as large as the maximum distance between two 
 * positions up to a specified cap
 */
class MinimumDistanceIndex {


    public: 
    ///Constructor for the distance index.
    ///Cap is the distance up to which the maximum distance will give a 
    ///reliable bound - if there is a path with length greater than cap, 
    ///then the maximum distance will be at least cap 
    ///If the cap is set to 0 (default), then the maximum distance index is not
    ///included
    MinimumDistanceIndex (const HandleGraph* graph, const SnarlManager* snarl_manager,
                            int64_t cap = 0);

    
    //Constructor to load index from serialization
    MinimumDistanceIndex (istream& in);
    
    //Default constructor; load() must be called next.
    MinimumDistanceIndex ();

    //Serialize object into out
    void serialize(ostream& out) const;

    //Load serialized object from in. Does not rely on the internal graph or 
    //snarl manager pointers.
    void load(istream& in);
    
    //Get the length of the given node
    int64_t node_length(id_t id) const;

    ///Get the minimum distance between two positions
    /// Distance includes only one of the positions. The distance from a 
    /// position to itself would be 1
    ///If there is no path between the two positions then the distance is -1
    int64_t min_distance( pos_t pos1, pos_t pos2) const;

    ///Get a maximum distance bound between the positions, ignoring direction
    ///Returns a positive value even if the two nodes are unreachable
    int64_t max_distance(pos_t pos1, pos_t pos2) const;


    //Given an alignment to a graph and a range, find the set of nodes in the
    //graph for which the minimum distance from the position to any position
    //in the node is within the given distance range
    //If look_forward is true, then start from the start of the path forward,
    //otherwise start from the end going backward
    void subgraph_in_range(const Path& path, const HandleGraph* super_graph, int64_t min_distance, int64_t max_distance, 
                         std::unordered_set<id_t>& sub_graph, bool look_forward);


    //Given a position, return distances that can be stored by a minimizer
    //
    //If the position is on a boundary node of a top level chain, then return true, and 
    //a unique identifier for the connected component that the node is on and
    //the offset of the position in the root chain - the minimum distance from the beginning of the chain to 
    //the position
    //The second bool will be false and the remaining size_t's will be 0
    //
    //If the position is on a child node of a top-level simple bubble (bubble has no children and nodes connect only to boundaries)
    //return false, 0, 0, true, and the rank of the bubble in its chain, the length of the start
    //node of the snarl, the length of the end node (relative to a fd traversal of the chain), and
    //the length of the node
    //
    //If the position is not on a root node (that is, a boundary node of a snarl in a root chain), returns
    //false and MIPayload::NO_VALUE for all values
    tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool>  get_minimizer_distances (pos_t pos);

    //What is the length of the top level chain that this node belongs to?
    int64_t top_level_chain_length(id_t node_id);

    size_t get_connected_component(id_t node_id);

    ///Helper function to find the minimum value that is not -1
    static int64_t min_pos(vector<int64_t> vals);

    ///Helper function to find the minimum value that is not -1
    static int64_t min_pos(int64_t x, int64_t y) {
        return static_cast<int64_t>(std::min(static_cast<uint64_t>(x), static_cast<uint64_t>(y)));
    }

    ///print the distance index for debugging
    void print_self();
    // Prints the number of nodes in each snarl netgraph and number of snarls in each chain
    void print_snarl_stats();

    protected:

    /** Index for calculating minimum distances among nodes in a snarl
     * Stores minimum distances between nodes in the netgraph of a snarl
     * Also keeps track of the parent of the snarl
    */
    class SnarlIndex {
        

        public:
        
            ///Constructor for the snarl index
            ///if inChain is true, parent_id and rev_in_parent are for the chain
            /// chain the snarl participates in. otherwise, for the parent snarl
            SnarlIndex(id_t parent_id, bool rev_in_parent, 
                       id_t id_in_parent, id_t end_id, bool is_unary_snarl, size_t depth,
                       size_t num_nodes, bool in_chain);
           
            //Construct an empty SnarlIndex. Must call load after construction to populate it 
            SnarlIndex();

            //Load data from serialization
            void load(istream& in, bool include_component);

            ///Serialize the snarl
            void serialize(ostream& out) const;
            
            ///Distance between start and end, not including the lengths of
            ///the two nodes
            ///start and end are the ranks of the node+direction, given by
            ///primary_snarls and secondary_snarl
            ///Only works for nodes heading their chains (which represent the 
            ///chains), or snarl boundaries.
            ///Rank 0 is the start node and rank num_nodes*2-1 is the end node,
            /// both pointing into the snarl
            int64_t snarl_distance(size_t start, size_t end) const {
                return int64_t(distances[index(start, end)]) - 1;
            }
             
            ///Length of a node in the netgraph of the snarl
            ///If it is a node, then the length of the node. If it is a snarl or
            ///chain, then the shortest distance between the boundaries
            /// i is the rank of the node in the snarl
            int64_t node_length(size_t i) const {
                return distances[i / 2] - 1;
            }
        
            ///Total length of the snarl-shortest distance from start to end
            ///including the lengths of boundary nodes
            int64_t snarl_length() const;

            ///Given distances from a position to either end of a node, find the
            ///shortest distance from that position to the start and end 
            ///nodes of the snarl
            ///rank is in the forward direction, but checks both forward and reverse
            
            pair<int64_t, int64_t> dist_to_ends(size_t rank, 
                                              int64_t distL, int64_t distR) const;

            ///For use during construction,
            ///add the distance from start to end to the index
            void insert_distance(size_t start, size_t end, int64_t dist);


            void print_self();

        protected:
 
            /// Store the distance between every pair nodes, not including the 
            /// lengths of the nodes. 
            /// The lengths of each of the nodes are stored as the first n entries
            /// Distances stored are 1 greater than actual distances 
            /// -1 (stored as 0) indicates no path
            /// For child snarls that are unary or only connected to one node
            /// in the snarl, distances between that node leaving the snarl
            /// and any other node is -1
            int_vector<> distances;

            ///True if this snarl is in a chain
            bool in_chain;

            /// id of the parent snarl or chain of this snarl 
            ///0 if this is a top level snarl with no chain
            id_t parent_id;

            bool rev_in_parent;

            ///id of this snarl in the parent. If the parent is a chain, then 
            ///the id of the boundary node that occurs first in the order of 
            ///the chain
            id_t id_in_parent;
            //id of the boundary node opposite id_in_parent
            id_t end_id;

            ///Number of nodes in the snarl
            size_t num_nodes;

            ///Depth in the snarl tree - 0 for root
            size_t depth;
            
            ///True if this snarl is a unary snarl
            ///Since the start and end node are the same, the last ranking
            ///node is no longer the end node if this is true
            bool is_unary_snarl;

            ///True if all children are nodes (not snarls or chains) and for every child node,
            //there are only edges to the boundary nodes
            bool is_simple_snarl;


            ///The maximum width of the snarl - the maximum of all minimum distance paths from each node to 
            //both ends of the snarl
            int64_t max_width;

            ///The index into distances for distance start->end
            size_t index(size_t start, size_t end) const;


        friend class MinimumDistanceIndex;
        friend class SnarlSeedClusterer;
    };


    /**Stores distances between snarls in a chain*/
    class ChainIndex {

        public:
        
            ///Constructor for ChainIndex
            /// loops is true if the chain loops - the start and end node are
            /// the same
            /// length is the number of snarls in the chain
            ChainIndex(id_t parent_id, id_t id_in_parent, id_t end_id, bool rev_in_parent,
                       bool loops, size_t length);

            //Constructor from vector of ints after serialization
            ChainIndex();
            //Load data from serialization
            void load(istream& in);

            ///Serialize the chain
            void serialize(ostream& out) const;
       
             
            ///Distance between two node sides in a chain. 
            ///size_t values specify the nodes - rank of the node in the chain, 
            ///bool specifies the direction the node is traversed i.e. if both
            /// bools are false, then the distance is between the left side of
            /// the start node traversing it forward to the left side of the 
            /// end node traversing forward. Orientation is relative to the 
            /// direction the chain is traversed in
            
            int64_t chain_distance(pair<size_t, bool> start, 
                                  pair<size_t, bool> end, int64_t startLen, 
                                  int64_t endLen, bool check_loop=false) const;


            //Length of entire chain
            int64_t chain_length() const {
                return prefix_sum[prefix_sum.size() - 1] - 1;
            }

            void print_self();

        protected:


            ///Dist from start of chain to start of each boundary node of
            ///all snarls in the chain
            ///The first value should be 0 according to this scheme but it is
            ///the length of the first node in the chain. Similarly, an extra
            ///value is stored at the end of the vector that is the length of the
            ///entire chain
            int_vector<> prefix_sum;

            ///For each boundary node of snarls in the chain, the distance
            /// from the start of the node traversing forward to the end of 
            /// the same node traversing backwards -directions relative to the 
            /// direction the node is traversed in the chain
            int_vector<> loop_fd;
    
            ///For each boundary node of snarls in the chain, the distance
            /// from the end of the node traversing backward to the start of 
            /// the same node traversing forward
            int_vector<> loop_rev;

            /// id of parent snarl of the chain 
            ///0 if top level chain
            id_t parent_id;

            bool rev_in_parent;

            //Id of the start node of this chain 
            id_t id_in_parent; 
            //id of the end node in the chain
            id_t end_id;

            //True if the chain loops - if the start and end node are the same
            bool is_looping_chain; 
            //Sum of all max widths of the snarls
            int64_t max_width;

            /// Helper function for chainDistance. Used to find the distance
            /// in a looping chain by taking the extra loop
            int64_t loop_distance(pair<size_t, bool> start, 
                                  pair<size_t, bool> end, int64_t startLen, 
                                  int64_t endLen) const;

        friend class MinimumDistanceIndex;   
        friend class SnarlSeedClusterer;
    }; 





    ///////// Data members of overall distance index


    id_t min_node_id; //minimum node id of the graph
    id_t max_node_id; //maximum node id of the graph

    ///vector of all SnarlIndex objects
    vector<SnarlIndex> snarl_indexes;

    ///vector of all ChainIndex objects
    vector< ChainIndex> chain_indexes;

    //Each connected component of the graph gets a unique identifier
    //Identifiers start at 1, 0 indicates that it is not in a component
    //Assigns each node to its connected component
    sdsl::int_vector<> node_to_component;
    //TODO: These could be one vector but they're small enough it probably doesn't matter
    sdsl::int_vector<> component_to_chain_length;
    sdsl::int_vector<> component_to_chain_index;

    //Keep track of the snarls ans chains each node belongs to
    enum SnarlType {TRIVIAL, SIMPLE, NONSIMPLE, NONE};
    //TODO: Might be easier to make things functions of the overall distance index- node length at least would make things simpler

    sdsl::int_vector<> node_assignments;



    ///The total depth of the snarl tree, starting from 0
    size_t tree_depth;


    ///True if we are including the maximum distance index
    bool include_maximum;

    //////Indexes for maximum distances
 
    ///For each node in the graph, store the minimum and maximum
    ///distances from a tip to the node
    sdsl::int_vector<> min_distances;
    sdsl::int_vector<> max_distances;


    //Header for the serialized file
    string file_header = "distance index version 2.2";
    //TODO: version 2 (no .anything) doesn't include component but we'll still accept it
    //version 2.1 doesn't include snarl index.is_simple_snarl and will break if we try to load it 
    bool include_component; //TODO: This is true for version 2.2 so it includes node_to_component, etc. 

    ////// Private helper functions
 




    ///Helper function for constructor - populate the minimum distance index
    ///Given the top level snarls
    //Returns the length of the chain
    int64_t calculate_min_index(const HandleGraph* graph, 
                      const SnarlManager* snarl_manager, const Chain* chain, 
                       size_t parent_id, bool rev_in_parent, 
                       bool trivial_chain, size_t depth, size_t component_num); 

    void populate_snarl_index(const HandleGraph* graph, const SnarlManager* snarl_manager, const NetGraph& ng,
                            const Snarl* snarl, bool snarl_rev_in_chain, size_t snarl_assignment, 
                            hash_set<pair<id_t, bool>>& all_nodes, size_t depth, size_t component_num);

    ///Compute min_distances and max_distances, which store
    /// distances needed for maximum distance calculation
    ///Only used if include_maximum is true
    void calculate_max_index(const HandleGraph* graph, int64_t cap); 



    ///Helper for subgraphInRange
    ///Given starting handles in the super graph and the distances to each handle (including the start position and 
    //the first position in the handle), add all nodes within the givendistance range to the subgraph
    //Ignore all nodes in seen_nodes (nodes that are too close)
    void add_nodes_in_range(const HandleGraph* super_graph, int64_t min_distance, int64_t max_distance, 
                         std::unordered_set<id_t>& sub_graph, vector<tuple<handle_t, int64_t>>& start_nodes,
                         hash_set<pair<id_t, bool>>& seen_nodes);

    ///Helper function for distance calculation
    ///Returns the distance to the start of and end of a node/snarl in
    //commonAncestor and the node id of the node in the common ancestor.
    ///The node in the common ancestor is a pair of <node id, rev>
    /// commonAncestor is a pair of the node_id and true if it is a chain
    ///rev is false if the pos is the start pos ( if it must be traversed 
    ///forward) and false if it is the end pos (if it must be reached in 
    ///the direction of pos)
    
    tuple<int64_t, int64_t, pair<id_t, bool>> dist_to_common_ancestor(
                pair<size_t, bool> common_ancestor, pos_t& pos, bool rev) const;



    int64_t get_primary_snarl_assignment  (id_t id) const { return (int64_t) node_assignments[(id-min_node_id) * 6 + 0] - 1;}
    int64_t get_primary_snarl_rank        (id_t id) const { return (int64_t) node_assignments[(id-min_node_id) * 6 + 1] - 1;}
    int64_t get_secondary_snarl_assignment(id_t id) const { return (int64_t) node_assignments[(id-min_node_id) * 6 + 2] - 1;}
    int64_t get_secondary_snarl_rank      (id_t id) const { return (int64_t) node_assignments[(id-min_node_id) * 6 + 3] - 1;}
    int64_t get_chain_assignment          (id_t id) const { return (int64_t) node_assignments[(id-min_node_id) * 6 + 4] - 1;}
    int64_t get_chain_rank                (id_t id) const { return (int64_t) node_assignments[(id-min_node_id) * 6 + 5] - 1;}

    void set_primary_snarl_assignment  (id_t id, int64_t val)    {node_assignments[(id-min_node_id) * 6 + 0] = val + 1;}
    void set_primary_snarl_rank        (id_t id, int64_t val)    {node_assignments[(id-min_node_id) * 6 + 1] = val + 1;}
    void set_secondary_snarl_assignment(id_t id, int64_t val)    {node_assignments[(id-min_node_id) * 6 + 2] = val + 1;}
    void set_secondary_snarl_rank      (id_t id, int64_t val)    {node_assignments[(id-min_node_id) * 6 + 3] = val + 1;}
    void set_chain_assignment          (id_t id, int64_t val)    {node_assignments[(id-min_node_id) * 6 + 4] = val + 1;}
    void set_chain_rank                (id_t id, int64_t val)    {node_assignments[(id-min_node_id) * 6 + 5] = val + 1;}

    friend class SnarlIndex;
    friend class ChainIndex;
    friend class SnarlSeedClusterer;


};

/**
 * The encoding of distances for positions in top-level chains or top-level simple bubbles.
 * Either stores (chain id, chain offset) for a position on a top-level chain, or
 * (snarl rank, node length, start length, end length) for a position on a simple bubble
 * We store this information in the minimizer index.
 */
/*
Simple bubble: 
    
 8 bit  |     1    |        24           |    10     |     10   |    10     |    1
  ---   |  is rev  | snarl rank in chain | start len | end len  | node len  |  is_node
   
Top level chain 
     
    31 bit   |    32    |     1
component id |  offset  |  is_node


is_node is true if it is a top-level chain node, false if it is a simple bubble
*/

struct MIPayload {
    typedef std::uint64_t code_type; // We assume that this fits into gbwtgraph::payload_type.

    constexpr static code_type NO_CODE = std::numeric_limits<code_type>::max();
    constexpr static size_t NO_VALUE = std::numeric_limits<size_t>::max(); // From offset_in_root_chain().

    constexpr static size_t NODE_LEN_OFFSET = 1;
    constexpr static size_t END_LEN_OFFSET = 11;
    constexpr static size_t START_LEN_OFFSET = 21;
    constexpr static size_t RANK_OFFSET = 31;
    constexpr static size_t REV_OFFSET = 55;

    
    constexpr static size_t LENGTH_WIDTH = 10;
    constexpr static size_t RANK_WIDTH = 24;
    constexpr static code_type LENGTH_MASK = (static_cast<code_type>(1) << LENGTH_WIDTH) - 1;
    constexpr static code_type RANK_MASK = (static_cast<code_type>(1) << RANK_WIDTH) - 1;
    

    
    constexpr static size_t ID_OFFSET = 33;
    constexpr static size_t ID_WIDTH = 31;
    constexpr static size_t OFFSET_WIDTH = 32;
    constexpr static code_type OFFSET_MASK = (static_cast<code_type>(1) << OFFSET_WIDTH) - 1;

    static code_type encode(std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_pos) {
        bool is_top_level_node; size_t component; size_t offset; //Values for a top level chain
        bool is_top_level_snarl; size_t snarl_rank; size_t node_length; size_t start_length; size_t end_length; bool is_rev;//values for a bubble
        std::tie(is_top_level_node, component, offset, is_top_level_snarl, snarl_rank, start_length, end_length, node_length, is_rev) = chain_pos;

        if (!is_top_level_node && ! is_top_level_snarl) {

            return NO_CODE;

        } else if (is_top_level_node) {
            //Top level node in chain

            if (component >= (static_cast<code_type>(1) << 31) - 1 
                || offset >= static_cast<size_t>(OFFSET_MASK) ) {
                //If the values are too large to be stored
                return NO_CODE;
            }

            return (component << ID_OFFSET) | (offset << 1) | static_cast<code_type>(true);

        } else {
            //Top level simple bubble

            if (snarl_rank >= static_cast<size_t>(RANK_MASK)
                || start_length >= static_cast<size_t>(LENGTH_MASK)
                || end_length >=  static_cast<size_t>(LENGTH_MASK)
                || node_length >= static_cast<size_t>(LENGTH_MASK) ){
                //If the values are too large to be stored
                return NO_CODE;
            }

            return (static_cast<code_type>(is_rev) << REV_OFFSET) | (snarl_rank << RANK_OFFSET) | (start_length << START_LEN_OFFSET) | (end_length << END_LEN_OFFSET) | (node_length << NODE_LEN_OFFSET) ;
        }
    }

    static std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> decode(code_type code) {
        if (code == NO_CODE) {
            return std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool>(false, NO_VALUE, NO_VALUE, false, NO_VALUE, NO_VALUE, NO_VALUE, NO_VALUE, false);
        } else if ((code & (static_cast<code_type>(1))) == (static_cast<code_type>(1))) {
            //This is a top-level chain
            return std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool>
                    (true, code >> ID_OFFSET, code >> 1 & OFFSET_MASK, false, NO_VALUE, NO_VALUE, NO_VALUE, NO_VALUE, false);
        } else {
            //This is a top-level bubble
            return std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool>
                    (false, NO_VALUE, NO_VALUE, true,
                      code >> RANK_OFFSET & RANK_MASK, 
                      code >> START_LEN_OFFSET & LENGTH_MASK, 
                      code >> END_LEN_OFFSET & LENGTH_MASK, 
                      code >> NODE_LEN_OFFSET & LENGTH_MASK,
                      code >> REV_OFFSET & static_cast<code_type>(1));
        }
    }
};

}

#endif
