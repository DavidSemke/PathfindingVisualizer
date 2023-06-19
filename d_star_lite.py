import min_heap
from search_tools import *

K_VARIABLE = 0

# k_variable is the variable described in the D* Lite algorithm which adjusts the heuristic key values
# in order to avoid reordering heap
def calculate_keys_d_lite(node, start, g_dict, rhs_dict):
    global ACCESSES

    g_value = g_dict[node]
    rhs_value = rhs_dict[node]

    key1 = (
        min(g_value, rhs_value) 
        + heuristic(start.get_pos(), node.get_pos()) 
        + K_VARIABLE
    )
    key2 = min(g_dict[node], rhs_dict[node])

    ACCESSES += 2

    return key1, key2


def update_node_d_lite(
        node_to_update, start, g_dict, rhs_dict, open_set
):
    global PERCOLATES
    global ACCESSES
    global COUNT

    index_of_item = None
    index = 0

    for item in open_set:
        
        if node_to_update is item[2]:
            index_of_item = index
            break
        
        index += 1

    locally_inconsistent = (
        g_dict[node_to_update] != rhs_dict[node_to_update]
    )

    ACCESSES += 2

    if locally_inconsistent and index_of_item is not None:
        k1, k2 = calculate_keys_d_lite(
            node_to_update, start, g_dict, rhs_dict
        )
        PERCOLATES += min_heap.heapremove(open_set, index_of_item)
        PERCOLATES += min_heap.heappush(
            open_set, ((k1, k2), COUNT, node_to_update)
        )
        COUNT += 1

    elif locally_inconsistent and index_of_item is None:
        k1, k2 = calculate_keys_d_lite(
            node_to_update, start, g_dict, rhs_dict
        )
        PERCOLATES += min_heap.heappush(
            open_set, ((k1, k2), COUNT, node_to_update)
        )
        COUNT += 1

        if not (
            node_to_update.is_invis_barrier() 
            or node_to_update.is_path()
        ):
            node_to_update.make_open()

    elif not locally_inconsistent and index_of_item is not None:
        PERCOLATES += min_heap.heapremove(open_set, index_of_item)


def d_star_lite_compute_shortest_path(
        draw_func, g_dict, rhs_dict, open_set, env
):
    global PERCOLATES
    global EXPANSIONS
    global ACCESSES
    global COUNT

    _, start, end = env

    while (
        (open_set and open_set[0][0] < calculate_keys_d_lite(
        start, start, g_dict, rhs_dict))
        or rhs_dict[start] != g_dict[start]
    ):
        current = open_set[0][2]
        k_old = open_set[0][0]
        k_new = calculate_keys_d_lite(current, start, g_dict, rhs_dict)

        ACCESSES += 2  # referring to those in while loop check

        if k_old < k_new:
            PERCOLATES += min_heap.heapreplace(
                open_set, (k_new, COUNT, current)
            )[1]
            COUNT += 1
        
        else:
            ACCESSES += 2
            
            if g_dict[current] > rhs_dict[current]:
                ACCESSES += 2
                g_dict[current] = rhs_dict[current]
                PERCOLATES += min_heap.heappop(open_set)[1]
                
                for node in current.neighbors:
                    
                    if node is not end:
                        rhs_dict[node] = min(
                            rhs_dict[node], 
                            g_dict[current] + NODE_TO_NODE_DIST
                        )
                        ACCESSES += 3

                    update_node_d_lite(
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
                        rhs_dict[node] == old_g_value + NODE_TO_NODE_DIST
                        and node is not end
                    ):
                        min_dist = float("inf")

                        for n in node.neighbors:
                            ACCESSES += 1
                            poss_new_min_dist = g_dict[n]
                            
                            if poss_new_min_dist < min_dist:
                                min_dist = poss_new_min_dist

                        ACCESSES += 1
                        rhs_dict[node] = min_dist + NODE_TO_NODE_DIST

                    update_node_d_lite(
                        node, start, g_dict, rhs_dict, open_set
                    )

        EXPANSIONS += 1

        draw_func()

    ACCESSES += 1

    return g_dict[start] != float('inf')


def perform_d_star_lite(draw_func, env):
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

    grid, start, end = env

    # priority queue as a heap
    open_set_heap = []
    # create dictionary of g_values set to infinity
    g_dict = {node: float("inf") for row in grid for node in row}
    # create dictionary of rhs_values set to infinity
    rhs_dict = {node: float("inf") for row in grid for node in row}
    # set rhs_value of goal node to 0
    rhs_dict[end] = 0
    ACCESSES += 1

    open_set_heap.append((
        (heuristic(end.get_pos(), start.get_pos()), 0), COUNT, end
    ))
    COUNT += 1

    origin = start
    true_origin = origin  # origin will be updated, so this is to remember first origin position

    # make all invis barriers within range of agent's vision at start position visible
    for n in true_origin.neighbors:
        
        if n.is_invis_barrier():
            n.make_vis_barrier()
            
            for n_neighbor in n.neighbors:
                n_neighbor.update_neighbors(grid)

    draw_func()

    # compute_shortest_path() returns True if a path, ultimately the shortest, to the goal is found
    if d_star_lite_compute_shortest_path(
        draw_func, g_dict, rhs_dict, open_set_heap, env
    ):
        
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

            # the next step simulates scanning for changes in 
            # edge costs
            # changes to graph can occur one node away from the 
            # start node in any direction
            nodes_changed = []
            
            for n in start.neighbors:
                
                if n.is_invis_barrier():
                    n.make_vis_barrier()
                    nodes_changed.append(n)
                    # remove from heap if present
                    index = 0
                    
                    for item in open_set_heap:
                        
                        if n is item[2]:
                            PERCOLATES += min_heap.heapremove(
                                open_set_heap, index
                            )
                            break
                        
                        index += 1

            if nodes_changed:
                K_VARIABLE += heuristic(
                    origin.get_pos(), start.get_pos()
                )
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
                            NODE_TO_NODE_DIST
                            and n_neighbor is not end
                        ):
                            min_dist = float("inf")
                            
                            for neighbor_of_n_neighbor in n_neighbor.neighbors:
                                ACCESSES += 1
                                poss_new_min_dist = g_dict[neighbor_of_n_neighbor]
                                
                                if poss_new_min_dist < min_dist:
                                    min_dist = poss_new_min_dist

                            ACCESSES += 1
                            rhs_dict[n_neighbor] = (
                                min_dist + NODE_TO_NODE_DIST
                            )

                        update_node_d_lite(
                            n_neighbor, start, g_dict, rhs_dict, 
                            open_set_heap
                        )

                    for node in [
                        node for row in grid for node in row
                    ]:
                        if node.is_closed() or node.is_open():
                            node.reset()

                if not d_star_lite_compute_shortest_path(
                    draw_func, g_dict, rhs_dict, open_set_heap,
                    env
                ):
                    # leave while loop if a path to start from goal 
                    # does not exist
                    break  

        if start is end:
            true_origin.make_original_start()
            print("Journey completed via D* Lite.")
        
        else:
            print("Journey unsuccessful - D* Lite failed to find path to goal.")
    else:
        print("Journey unsuccessful - D* Lite failed to find path to goal.")

    print("Heap percolates: " + str(PERCOLATES))
    print("Node expansions: " + str(EXPANSIONS))
    print("Node accesses: " + str(ACCESSES))