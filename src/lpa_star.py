import min_heap as heap
from grid import reset_nodes
from search_tools import *

# nodes are prioritized by their keys, where smaller key values are 
# of greater priority
def lpa_calculate_keys(node, end, g_dict, rhs_dict):
    global ACCESSES

    g_value = g_dict[node]
    rhs_value = rhs_dict[node]
    
    # f_value correspondence
    h = heuristic(end.get_pos(), node.get_pos())
    key1 = min(g_value, rhs_value) + h
    
    # g_value correspondence
    key2 = min(g_value, rhs_value)

    ACCESSES += 2

    return key1, key2


def lpa_update_node(node_to_update, end, g_dict, rhs_dict, open_set):
    global COUNT
    global PERCOLATES
    global ACCESSES

    item_index = None
    
    for i, item in enumerate(open_set):
        
        if node_to_update is item[2]:
            item_index = i
            break

    locally_inconsistent = (
        g_dict[node_to_update] != rhs_dict[node_to_update]
    )

    ACCESSES += 2

    if locally_inconsistent and item_index is not None:
        k1, k2 = lpa_calculate_keys(
            node_to_update, end, g_dict, rhs_dict
        )
        PERCOLATES += heap.heapremove(open_set, item_index)
        PERCOLATES += heap.heappush(
            open_set, ((k1, k2), COUNT, node_to_update)
        )
        COUNT += 1

    elif locally_inconsistent and item_index is None:
        k1, k2 = lpa_calculate_keys(
            node_to_update, end, g_dict, rhs_dict
        )
        PERCOLATES += heap.heappush(
            open_set, ((k1, k2), COUNT, node_to_update)
        )
        COUNT += 1

        if not node_to_update.is_invis_barrier():
            node_to_update.make_open()

    elif not (locally_inconsistent or item_index is None):
        PERCOLATES += heap.heapremove(open_set, item_index)


#  draw_func is necessary for updating grid
#  start and end are start and goal nodes
def lpa_shortest_path(
        draw_func, g_dict, rhs_dict, open_set, start, end
):
    global PERCOLATES
    global EXPANSIONS
    global ACCESSES

    # get top item in priority queue
    while (
        (open_set and open_set[0][0] < lpa_calculate_keys(
        end, end, g_dict, rhs_dict))
        or rhs_dict[end] != g_dict[end]
    ):
        current = open_set[0][2]

        # accesses in while loop and in following if statement
        ACCESSES += 4

        if g_dict[current] > rhs_dict[current]:
            g_dict[current] = rhs_dict[current]
            PERCOLATES += heap.heappop(open_set)[1]
            ACCESSES += 2
            
            for node in current.neighbors:
                
                if node is not start:
                    rhs_dict[node] = (
                        min(rhs_dict[node], 
                        g_dict[current] + EDGE_COST)
                    )
                    ACCESSES += 3

                lpa_update_node(
                    node, end, g_dict, rhs_dict, open_set
                )

            if (
                current is not start 
                and not current.is_invis_barrier()
            ):
                current.make_closed()

        else:
            old_g_value = g_dict[current]
            g_dict[current] = float('inf')
            ACCESSES += 2

            neighbors_and_current = current.neighbors + [current]
            
            for node in neighbors_and_current:
                ACCESSES += 1
                
                if (
                    (rhs_dict[node] == old_g_value +
                    EDGE_COST or node is current)
                    and node is not start
                ):
                    min_dist = float("inf")
                    
                    for n in node.neighbors:
                        ACCESSES += 1
                        poss_new_min_dist = g_dict[n]
                        
                        if poss_new_min_dist < min_dist:
                            min_dist = poss_new_min_dist

                    rhs_dict[node] = min_dist + EDGE_COST
                    ACCESSES += 1

                lpa_update_node(node, end, g_dict, rhs_dict, open_set)

        EXPANSIONS += 1

        draw_func()

    ACCESSES += 1

    if g_dict[end] == float('inf'): return False
    
    reconstruct_path(start, end, g_dict, draw_func)
    end.make_end()
    start.make_start()
         
    return True


# main method for LPA*
def lpa_star(draw_func, env):
    global PERCOLATES
    global EXPANSIONS
    global ACCESSES
    global COUNT
    PERCOLATES = 0
    EXPANSIONS = 0
    ACCESSES = 0
    COUNT = 0

    grid = env['grid']
    start = env['start']
    end = env['end']
    invis_barriers = env['invis_barriers']

    open_set_heap = []
    # create dictionary of g_values set to infinity
    g_dict = {
        node: float("inf") for row in grid for node in row
    }
    # create dictionary of rhs_values set to infinity
    rhs_dict = {
        node: float("inf") for row in grid for node in row
    }

    # set rhs_value of start node to 0
    rhs_dict[start] = 0
    ACCESSES += 1

    h = heuristic(end.get_pos(), start.get_pos())
    open_set_heap.append(((h, 0), COUNT, start))
    COUNT += 1

    lpa_shortest_path(
        draw_func, g_dict, rhs_dict, open_set_heap, start, end
    )

    for b in invis_barriers:
        b.make_vis_barrier()
        
        for i, item in enumerate(open_set_heap):
            
            if not b is item[2]: continue
            
            PERCOLATES += heap.heapremove(
                open_set_heap, i
            )
            
            break

        b.update_neighbors(grid)

        for n_neighbor in b.neighbors:
            n_neighbor.update_neighbors(grid)
            ACCESSES += 2
            
            if (
                rhs_dict[n_neighbor] == g_dict[b] + EDGE_COST
                and n_neighbor is not start
            ):
                min_dist = float("inf")
                
                for neighbor_of_n_neighbor in n_neighbor.neighbors:
                    ACCESSES += 1
                    poss_new_min_dist = g_dict[neighbor_of_n_neighbor]
                    
                    if poss_new_min_dist < min_dist:
                        min_dist = poss_new_min_dist
                
                rhs_dict[n_neighbor] = min_dist + EDGE_COST
                ACCESSES += 1

            lpa_update_node(
                n_neighbor, end, g_dict, rhs_dict, open_set_heap
            )

        # clear previous activity
        l = lambda n: n.is_path() or n.is_closed() or n.is_open()
        reset_nodes(grid, l)
            
        draw_func()

        lpa_shortest_path(
            draw_func, g_dict, rhs_dict, open_set_heap, start, end
        )

    print("LPA* completed trials.")
    
    print("Heap percolates: " + str(PERCOLATES))
    print("Node expansions: " + str(EXPANSIONS))
    print("Node accesses: " + str(ACCESSES))