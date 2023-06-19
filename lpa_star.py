import min_heap
from search_tools import *

# nodes are prioritized by their keys, where smaller key values are 
# of greater priority
def calculate_keys_lpa(node, end, g_dict, rhs_dict):
    global ACCESSES

    g_value = g_dict[node]
    rhs_value = rhs_dict[node]
    
    # f_value correspondence
    key1 = (
        min(g_value, rhs_value) 
        + heuristic(end.get_pos(), node.get_pos())
    )
    
    # g_value correspondence
    key2 = min(g_value, rhs_value)

    ACCESSES += 2

    return key1, key2


def update_node_lpa(node_to_update, end, g_dict, rhs_dict, open_set):
    global COUNT
    global PERCOLATES
    global ACCESSES

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
        k1, k2 = calculate_keys_lpa(
            node_to_update, end, g_dict, rhs_dict
        )
        PERCOLATES += min_heap.heapremove(open_set, index_of_item)
        PERCOLATES += min_heap.heappush(
            open_set, ((k1, k2), COUNT, node_to_update)
        )
        COUNT += 1

    elif locally_inconsistent and index_of_item is None:
        k1, k2 = calculate_keys_lpa(
            node_to_update, end, g_dict, rhs_dict
        )
        PERCOLATES += min_heap.heappush(
            open_set, ((k1, k2), COUNT, node_to_update)
        )
        COUNT += 1

        if not node_to_update.is_invis_barrier():
            node_to_update.make_open()

    elif not locally_inconsistent and index_of_item is not None:
        PERCOLATES += min_heap.heapremove(open_set, index_of_item)


#  draw_func is necessary for updating grid
#  start and end are start and goal nodes
def lpa_star_compute_shortest_path(
        draw_func, g_dict, rhs_dict, open_set, env
):
    global PERCOLATES
    global EXPANSIONS
    global ACCESSES

    _, start, end = env

    # get top item in priority queue
    while (
        (open_set and open_set[0][0] < calculate_keys_lpa(
        end, end, g_dict, rhs_dict))
        or rhs_dict[end] != g_dict[end]
    ):
        current = open_set[0][2]

        # accesses in while loop and in following if statement
        ACCESSES += 4

        if g_dict[current] > rhs_dict[current]:
            g_dict[current] = rhs_dict[current]
            PERCOLATES += min_heap.heappop(open_set)[1]
            ACCESSES += 2
            
            for node in current.neighbors:
                
                if node is not start:
                    rhs_dict[node] = (
                        min(rhs_dict[node], 
                        g_dict[current] + NODE_TO_NODE_DIST)
                    )
                    ACCESSES += 3

                update_node_lpa(
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
                    NODE_TO_NODE_DIST or node is current)
                    and node is not start
                ):
                    min_distance = float("inf")
                    
                    for n in node.neighbors:
                        ACCESSES += 1
                        possible_new_min_distance = g_dict[n]
                        
                        if possible_new_min_distance < min_distance:
                            min_distance = possible_new_min_distance

                    rhs_dict[node] = min_distance + NODE_TO_NODE_DIST
                    ACCESSES += 1

                update_node_lpa(node, end, g_dict, rhs_dict, open_set)

        EXPANSIONS += 1

        draw_func()

    ACCESSES += 1

    if g_dict[end] == float('inf'): return False
    
    reconstruct_path(env, g_dict, draw_func)
    end.make_end()
    start.make_start()
         
    return True


# main method for LPA*
def perform_lpa_star(draw_func, env, invis_barriers):
    global PERCOLATES
    global EXPANSIONS
    global ACCESSES
    global COUNT
    PERCOLATES = 0
    EXPANSIONS = 0
    ACCESSES = 0
    COUNT = 0

    grid, start, end = env

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

    open_set_heap.append((
        (heuristic(end.get_pos(), start.get_pos()), 0), COUNT, start
    ))
    COUNT += 1

    invis_barriers_index = 0

    lpa_star_compute_shortest_path(
        draw_func, g_dict, rhs_dict, open_set_heap, env
    )

    while invis_barriers_index < len(invis_barriers):
        b = invis_barriers[invis_barriers_index]
        b.make_vis_barrier()
        invis_barriers_index += 1
        index = 0
        
        for item in open_set_heap:
            
            if b is item[2]:
                PERCOLATES += min_heap.heapremove(
                    open_set_heap, index
                )
                break
            
            index += 1

        b.update_neighbors(grid)

        for n_neighbor in b.neighbors:
            n_neighbor.update_neighbors(grid)
            ACCESSES += 2
            
            if (
                rhs_dict[n_neighbor] == g_dict[b] + NODE_TO_NODE_DIST
                and n_neighbor is not start
            ):
                min_distance = float("inf")
                
                for neighbor_of_n_neighbor in n_neighbor.neighbors:
                    ACCESSES += 1
                    possible_new_min_distance = g_dict[neighbor_of_n_neighbor]
                    
                    if possible_new_min_distance < min_distance:
                        min_distance = possible_new_min_distance
                
                rhs_dict[n_neighbor] = min_distance + NODE_TO_NODE_DIST
                ACCESSES += 1

            update_node_lpa(
                n_neighbor, end, g_dict, rhs_dict, open_set_heap
            )

        # clear previous activity
        for node in [node for row in grid for node in row]:
            if node.is_path() or node.is_closed() or node.is_open():
                node.reset()

        draw_func()

        lpa_star_compute_shortest_path(
            draw_func, g_dict, rhs_dict, open_set_heap, env
        )

    else:
        print("LPA* completed trials.")

    print("Heap percolates: " + str(PERCOLATES))
    print("Node expansions: " + str(EXPANSIONS))
    print("Node accesses: " + str(ACCESSES))