import numpy as np


def find_contributing_cells(aCellID, aCellID_downslope, lCellID):
    """
    Find all cells that contribute (drain into) the given cell lCellID,
    including the cell itself.

    This is a BFS/DFS traversal of the upslope graph starting from lCellID.
    For repeated calls over all cells, prefer build_contributing_cells_all()
    which computes the full mapping in O(N) instead of O(N^2).

    Args:
        aCellID (array-like): Array of cell IDs.
        aCellID_downslope (array-like): Array of downstream cell IDs (same order as aCellID).
        lCellID: The target cell ID whose contributing cells are sought.

    Returns:
        tuple: (aCellID_contribution, aCellIndex_contribution)
            - aCellID_contribution: list of contributing cell IDs
            - aCellIndex_contribution: list of contributing cell indices
    """
    aCellID = np.asarray(aCellID)
    aCellID_downslope = np.asarray(aCellID_downslope)

    # Build upslope adjacency list once for this call
    id_to_index = {cid: idx for idx, cid in enumerate(aCellID)}

    # upslope_adj[i] = list of cell indices that drain directly into cell i
    nCell = len(aCellID)
    upslope_adj = [[] for _ in range(nCell)]
    for i in range(nCell):
        ds_id = aCellID_downslope[i]
        if ds_id in id_to_index:
            j = id_to_index[ds_id]
            upslope_adj[j].append(i)

    # BFS from lCellID upward through the upslope graph
    start_index = id_to_index[lCellID]
    visited = [start_index]
    queue = [start_index]
    head = 0
    while head < len(queue):
        current = queue[head]
        head += 1
        for upstream in upslope_adj[current]:
            visited.append(upstream)
            queue.append(upstream)

    aCellIndex_contribution = visited
    aCellID_contribution = [aCellID[idx] for idx in visited]

    return aCellID_contribution, aCellIndex_contribution


def build_contributing_cells_all(aCellID, aCellID_downslope):
    """
    Compute contributing cell index lists for ALL cells in O(N) using a
    single topological sort (Kahn's algorithm).

    This replaces calling find_contributing_cells() N times (which is O(N^2)).

    Args:
        aCellID (array-like): Array of cell IDs, length N.
        aCellID_downslope (array-like): Array of downstream cell IDs, length N.

    Returns:
        list[list[int]]: aCellIndex_contribution_all[i] is the list of cell
                         indices (0-based) that drain into cell i, including i itself.
    """
    aCellID = np.asarray(aCellID)
    aCellID_downslope = np.asarray(aCellID_downslope)
    nCell = len(aCellID)

    # Map cell ID -> index
    id_to_index = {cid: idx for idx, cid in enumerate(aCellID)}

    # Build upslope adjacency and compute in-degree (number of upstream neighbours)
    upslope_adj = [[] for _ in range(nCell)]
    in_degree = np.zeros(nCell, dtype=np.int64)

    for i in range(nCell):
        ds_id = aCellID_downslope[i]
        if ds_id in id_to_index:
            j = id_to_index[ds_id]
            upslope_adj[j].append(i)
            in_degree[j] += 1

    # Each cell starts by contributing to itself
    # Use lists for O(1) extend; convert to np.array at the end if needed
    contrib = [[] for _ in range(nCell)]
    for i in range(nCell):
        contrib[i].append(i)

    # Kahn's topological sort: process headwater cells first
    from collections import deque
    queue = deque(int(i) for i in np.where(in_degree == 0)[0])

    while queue:
        i = queue.popleft()
        ds_id = aCellID_downslope[i]
        if ds_id in id_to_index:
            j = id_to_index[ds_id]
            # Propagate i's full contributor list downstream to j
            contrib[j].extend(contrib[i])
            in_degree[j] -= 1
            if in_degree[j] == 0:
                queue.append(j)

    return contrib
