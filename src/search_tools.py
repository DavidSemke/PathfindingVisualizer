# universal edge cost (edges costing infinity are simply ignored)
EDGE_COST = 1

'''The following are global counters.'''

# PERCOLATES - Total heap percolates
# EXPANSIONS - Total expansions of nodes (expansion = visiting a node AND updating its g_value)
# ACCESSES - Total accesses of nodes (updating and reading node attributes, including g-values, f_values,
# and rhs-values)
# COUNT - Ensures each entry into the min heap is unique, so that a tie breaker between priorities leads to the oldest
# entry being prioritized

# Each of these variables reset to 0 at the start of an algorithm
PERCOLATES = 0
EXPANSIONS = 0
ACCESSES = 0
COUNT = 0


# Define the heuristic for a node. Uses Manhattan distance.
# Note that the order of points as parameters does not matter
def heuristic(p1, p2):  
    x1, y1 = p1
    x2, y2 = p2
    
    return abs(x1 - x2) + abs(y1 - y2)


# to reconstruct path, find the predecessor which minimizes g_value + edge cost, and
# repeat until the predecessor is the start node
# edge cost is fixed at NODE_TO_NODE_DISTANCE, so can just compare g-values
def reconstruct_path(start, end, g_dict, draw_func):
    current = end
    next_node = None

    while current is not start:
        min_dist = float("inf")

        for n in current.neighbors:
            poss_new_min_dist = g_dict[n]
            
            if poss_new_min_dist < min_dist:
                min_dist = poss_new_min_dist
                next_node = n
        
        current = next_node
        
        if not current.is_invis_barrier():
            current.make_path()
        
        draw_func()