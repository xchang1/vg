///
///  \file handle_snarls.hpp
///
///  Contains a handle-based interface for snarls.
///

#ifndef VG_HANDLE_SNARLS_HPP_INCLUDED
#define VG_HANDLE_SNARLS_HPP_INCLUDED

#include <iostream>
#include <functional>
#include "handle.hpp"

namespace vg {

/// Represents a snarl
struct snarl_handle_t { char data[2 * sizeof(nid_t)]; };
/// Represents a chain
struct chain_handle_t { char data[2 * sizeof(nid_t)]; };

/**
 * Represents a breaking-up of a graph into a bilayered tree of snarls and
 * chains, rooted at boundary-less root snarls.
 */
class SnarlChainDecomposition {
public:

    //////////
    // Tree traversal operations
    //////////
   
    /**
     * Loop over all root, boundary-less snarls in the decomposition.
     * One exists per connected component in the graph.
     */
    virtual void for_each_root_snarl(const function<bool(const snarl_handle_t&)>& iteratee) const = 0;
   
    /**
     * Loop over all the chains in a snarl, including single-node "terminal" or
     * "empty" chains.
     *
     * If the snarl is a root snarl, may include circular chains, which start
     * and end at the same node but contain child snarls.
     */
    virtual void for_each_child(const snarl_handle_t& snarl, const function<bool(const chain_handle_t&)>& iteratee) const = 0;
    
    /**
     * Loop over all the snarls in a chain, in order. All snarls are taken to
     * be forward in their partent chains.
     */
    virtual void for_each_child(const chain_handle_t& chain, const function<bool(const snarl_handle_t&)>& iteratee) const = 0;
    
    /**
     * Get the parent chain of any snarl that is not a root snarl.
     * May not be called for root snarls.
     */
    virtual chain_handle_t get_parent(const snarl_handle_t& snarl) const = 0;
    
    /**
     * Get the parent snarl of any chain.
     */
    virtual snarl_handle_t get_parent(const chain_handle_t& chain) const = 0;
    
    /**
     * Return true if there is another snarl after this one in its parent chain.
     */
    virtual bool has_next_in_chain(const snarl_handle_t& snarl) const = 0;
    
    /**
     * Return the next snarl in the same chain as the given snarl. Said next
     * snarl must exist.
     */
    virtual snarl_handle_t get_next_in_chain(const snarl_handle_t& snarl) const = 0;
    
    /**
     * Return true if there is another snarl before this one in its parent chain.
     */
    virtual bool has_prev_in_chain(const snarl_handle_t& snarl) const = 0;
    
    /**
     * Return the previous snarl in the same chain as the given snarl. Said
     * previous snarl must exist.
     */
    virtual snarl_handle_t get_prev_in_chain(const snarl_handle_t& snarl) const = 0;
    
    //////////
    // Snarl accessors
    //////////
    
    /**
     * Return true if a snarl is a boundary-less root snarl, and false
     * otherwise.
     */
    virtual bool is_root(const snarl_handle_t& snarl) const = 0;
    
    /**
     * Get the inward, start handle of a snarl, if it is not a root snarl.
     */
    virtual handle_t get_start(const snarl_handle_t& snarl) const = 0;
    
    /**
     * Get the outward, end handle of a snarl, if it is not a root snarl.
     */
    virtual handle_t get_end(const snarl_handle_t& snarl) const = 0;
    
    //////////
    // Chain accessors
    //////////
    
    /**
     * Return true if a chain is circular, and false otherwise.
     */
    virtual bool is_circular(const chain_handle_t& chain) const = 0;
    
    /**
     * Return true if a chain is a single node, and false otherwise.
     */
    virtual bool is_terminal(const chain_handle_t& chain) const = 0;
    
    /**
     * For a terminal chain, return the node it represents.
     */
    virtual handle_t get_terminal(const chain_handle_t& chain) const = 0;
    
    //////////
    // Index operations
    //////////
    
    /**
     * Find the lowest snarl containing the given handle. Will be the snarl the
     * handle reads into, if the handle's node is a snarl boundary, or the
     * snarl containing the terminal chain representing the handle, if the
     * handle's node is not a snarl boundary.
     */
    virtual snarl_handle_t find_containing_snarl(const handle_t& query) const = 0;
    
    /**
     * Find the lowest chain containing the given handle. Will be the chain
     * containing the snarl(s) that the handle's node bounds, if the handle's
     * node is a snarl boundary, or the terminal chain representing the handle,
     * otherwise.
     */
    virtual chain_handle_t find_containing_chain(const handle_t& query) const = 0;
   
    
    
    
};


}

#endif
