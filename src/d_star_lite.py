import min_heap as heap
from grid import reset_nodes
from search_tools import *

# k_variable is the variable described in the D* Lite algorithm which 
# adjusts the heuristic key values in order to avoid reordering heap
K_VARIABLE = 0

def d_lite_calculate_keys(node, start, g_dict, rhs_dict):
    global ACCESSES

    g_value = g_dict[node]
    rhs_value = rhs_dict[node]

    h = heuristic(start.get_pos(), node.get_pos())
    key1 = min(g_value, rhs_value) + h + K_VARIABLE
    key2 = min(g_dict[node], rhs_dict[node])

    ACCESSES += 2

    return key1, key2


def d_lite_update_node(
        node_to_update, start, g_dict, rhs_dict, open_set
):
    global PERCOLATES
    global ACCESSES
    global COUNT

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
        k1, k2 = d_lite_calculate_keys(
            node_to_update, start, g_dict, rhs_dict
        )
        PERCOLATES += heap.heapremove(open_set, item_index)
        PERCOLATES += heap.heappush(
            open_set, ((k1, k2), COUNT, node_to_update)
        )
        COUNT += 1

    elif locally_inconsistent and item_index is None:
        k1, k2 = d_lite_calculate_keys(
            node_to_update, start, g_dict, rhs_dict
        )
        PERCOLATES += heap.heappush(
            open_set, ((k1, k2), COUNT, node_to_update)
        )
        COUNT += 1

        if not (
            node_to_update.is_invis_barrier() 
            or node_to_update.is_path()
        ):
            node_to_update.make_open()

    elif not (locally_inconsistent or item_index is None):
        PERCOLATES += heap.heapremove(open_set, item_index)


def d_lite_shortest_path(
        draw_func, g_dict, rhs_dict, open_set, start, end
):
    global PERCOLATES
    global EXPANSIONS
    global ACCESSES
    global COUNT

    while (
        (open_set and open_set[0][0] < d_lite_calculate_keys(
        start, start, g_dict, rhs_dict))
        or rhs_dict[start] != g_dict[start]
    ):
        current = open_set[0][2]
        k_old = open_set[0][0]
        k_new = d_lite_calculate_keys(current, start, g_dict, rhs_dict)

        ACCESSES += 2  # referring to those in while loop check

        if k_old < k_new:
            PERCOLATES += heap.heapreplace(
                open_set, (k_new, COUNT, current)
            )[1]
            COUNT += 1
        
        else:
            ACCESSES += 2
            
            if g_dict[current] > rhs_dict[current]:
                ACCESSES += 2
                g_dict[current] = rhs_dict[current]
                PERCOLATES += heap.heappop(open_set)[1]
                
                for node in current.neighbors:
                    
                    if node is not end:
                        rhs_dict[node] = min(
                            rhs_dict[node], 
                            g_dict[current] + EDGE_COST
                        )
                        ACCESSES += 3

                    d_lite_update_node(
                        node, start, g_dict, rhs_dict, open_set
                    )

                if not (
                    current is end 
                    or current.is_invis_barrier() 
                    or current.is_path()
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
                        rhs_dict[node] == old_g_value + EDGE_COST
                        and node is not end
                    ):
                        min_dist = float("inf")

                        for n in node.neighbors:
                            ACCESSES += 1
                            poss_new_min_dist = g_dict[n]
                            
                            if poss_new_min_dist < min_dist:
                                min_dist = poss_new_min_dist

                        ACCESSES += 1
                        rhs_dict[node] = min_dist + EDGE_COST

                    d_lite_update_node(
                        node, start, g_dict, rhs_dict, open_set
                    )

        EXPANSIONS += 1

        draw_func()

    ACCESSES += 1

    return g_dict[start] != float('inf')


def d_star_lite(draw_func, env):
    global PERCOLATES
    global EXPANSIONS
    global ACCESSES
    global COUNT
    global K_VARIABLE
    PERCOLATES = 0
    EXPANSIONS = 0
    ACCESSES = 0
    COUNT = 0
    K_VARIABLE = 0

    grid = env['grid']
    start = env['start']
    end = env['end']

    # priority queue as a heap
    open_set_heap = []
    # create dictionary of g_values set to infinity
    g_dict = {node: float("inf") for row in grid for node in row}
    # create dictionary of rhs_values set to infinity
    rhs_dict = {node: float("inf") for row in grid for node in row}
    # set rhs_value of goal node to 0
    rhs_dict[end] = 0
    ACCESSES += 1

    h = heuristic(end.get_pos(), start.get_pos())
    open_set_heap.append(((h, 0), COUNT, end))
    COUNT += 1

    origin = start
    # origin will be updated, so this is to remember first origin position
    true_origin = origin  

    # make all invis barriers within range of agent's vision at start position visible
    for n in true_origin.neighbors:
        
        if not n.is_invis_barrier(): continue
        
        n.make_vis_barrier()
        
        for n_neighbor in n.neighbors:
            n_neighbor.update_neighbors(grid)

    draw_func()

    SUCCESS = "Journey completed via D* Lite."
    FAIL = "Journey unsuccessful - D* Lite failed to find path to goal."

    path_exists = d_lite_shortest_path(
        draw_func, g_dict, rhs_dict, open_set_heap, start, end
    )

    if not path_exists:
        print(FAIL)
        return

    while start is not end:
        next_start_node = None
        min_dist = float("inf")
        
        for n in start.neighbors:
            ACCESSES += 1
            poss_new_min_dist = g_dict[n]
            
            if poss_new_min_dist < min_dist:
                min_dist = poss_new_min_dist
                next_start_node = n

        start.make_path()
        start = next_start_node
        start.make_start()

        draw_func()

        # the next step simulates scanning for changes in edge costs
        # changes to graph can occur one node away from the start
        # node in any direction
        nodes_changed = []
        
        for n in start.neighbors:
            
            if not n.is_invis_barrier(): continue
            
            n.make_vis_barrier()
            nodes_changed.append(n)
            
            # remove from heap if present                    
            for i, item in enumerate(open_set_heap):
                
                if not n is item[2]: continue
                
                PERCOLATES += heap.heapremove(
                    open_set_heap, i
                )
                
                break

        if not nodes_changed: continue
        
        K_VARIABLE += heuristic(origin.get_pos(), start.get_pos())
        
        origin = start

        # note that some code has been omitted here, as the 
        # code would not apply to an environment with
        # solely traversable and not traversable edges (edge 
        # is either some constant or infinity)
        for n in nodes_changed:
            n.update_neighbors(grid)
            
            for n_neighbor in n.neighbors:
                n_neighbor.update_neighbors(grid)
                ACCESSES += 2
                
                if (
                    rhs_dict[n_neighbor] == g_dict[n] +
                    EDGE_COST
                    and n_neighbor is not end
                ):
                    min_dist = float("inf")
                    
                    for neighbor_of_n_neighbor in n_neighbor.neighbors:
                        ACCESSES += 1
                        poss_new_min_dist = g_dict[neighbor_of_n_neighbor]
                        
                        if poss_new_min_dist < min_dist:
                            min_dist = poss_new_min_dist

                    ACCESSES += 1
                    rhs_dict[n_neighbor] = min_dist + EDGE_COST

                d_lite_update_node(
                    n_neighbor, start, g_dict, rhs_dict, 
                    open_set_heap
                )

            l = lambda n: n.is_closed() or n.is_open()
            reset_nodes(grid, l)

        path_exists = d_lite_shortest_path(
            draw_func, g_dict, rhs_dict, open_set_heap,
            start, end
        )

        if not path_exists: break

    if start is end:
        true_origin.make_original_start()
        print(SUCCESS)
    
    else:
        print(FAIL)

    print("Heap percolates: " + str(PERCOLATES))
    print("Node expansions: " + str(EXPANSIONS))
    print("Node accesses: " + str(ACCESSES))