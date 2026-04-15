# Shortest Path Algorithms

**Computer Science Fundamentals Series**

Dijkstra · Bellman-Ford · Floyd-Warshall · A* search · Negative cycles · Routing algorithms

*Mid-level software engineer track -- 20 slides*

---

## Table of Contents

1. [Single-Source Shortest Path Problem](#slide-02--single-source-shortest-path-problem)
2. [BFS for Unweighted Graphs](#slide-03--bfs-for-unweighted-graphs)
3. [Dijkstra's Algorithm](#slide-04--dijkstras-algorithm)
4. [Dijkstra -- Priority Queue Implementation](#slide-05--dijkstra--priority-queue-implementation)
5. [Dijkstra -- Proof of Correctness](#slide-06--dijkstra--proof-of-correctness)
6. [Bellman-Ford Algorithm](#slide-07--bellman-ford-algorithm)
7. [Bellman-Ford -- Negative Cycle Detection](#slide-08--bellman-ford--negative-cycle-detection)
8. [SPFA -- Shortest Path Faster Algorithm](#slide-09--spfa--shortest-path-faster-algorithm)
9. [DAG Shortest Paths via Topological Sort](#slide-10--dag-shortest-paths-via-topological-sort)
10. [Floyd-Warshall -- All-Pairs Shortest Paths](#slide-11--floyd-warshall--all-pairs-shortest-paths)
11. [Johnson's Algorithm](#slide-12--johnsons-algorithm)
12. [A* Search](#slide-13--a-search)
13. [A* -- Heuristics & Admissibility](#slide-14--a--heuristics--admissibility)
14. [Bidirectional Search](#slide-15--bidirectional-search)
15. [Contraction Hierarchies](#slide-16--contraction-hierarchies)
16. [Applications -- Routing & Network Protocols](#slide-17--applications--routing--network-protocols)
17. [Applications -- Game AI & Beyond](#slide-18--applications--game-ai--beyond)
18. [Complexity Comparison](#slide-19--complexity-comparison)
19. [Summary & Further Reading](#slide-20--summary--further-reading)

---

## Slide 02 -- Single-Source Shortest Path Problem

### Definition

Given a weighted directed graph `G = (V, E)` with weight function `w: E -> R` and a source vertex `s`, find the shortest-weight path from `s` to every reachable vertex `v in V`.

- **Path weight** -- sum of edge weights along the path: `w(p) = sum w(v_i, v_{i+1})`
- **Shortest-path weight** -- `delta(s, v) = min { w(p) : s ~> v }` or `+inf` if unreachable
- **Predecessor sub-graph** -- encodes shortest-path tree rooted at `s`

### Variants

| Variant | Description |
|---------|------------|
| **Single-source** | One source to all vertices -- Dijkstra, Bellman-Ford |
| **Single-destination** | All vertices to one target -- reverse all edges, run single-source |
| **Single-pair** | One source to one target -- A*, bidirectional search |
| **All-pairs** | Every vertex to every vertex -- Floyd-Warshall, Johnson's |

### Optimal sub-structure

Any sub-path of a shortest path is itself a shortest path. This property is the foundation for every algorithm in this deck.

---

## Slide 03 -- BFS for Unweighted Graphs

BFS solves single-source shortest paths when every edge has weight 1. It processes vertices in order of increasing distance from the source.

```
BFS(G, s):
  for each vertex u in V:
    dist[u] = inf
  dist[s] = 0
  Q = empty queue
  enqueue(Q, s)
  while Q is not empty:
    u = dequeue(Q)
    for each edge (u, v) in E:
      if dist[v] == inf:
        dist[v] = dist[u] + 1
        prev[v] = u
        enqueue(Q, v)
```

### Properties

- **Time complexity** -- `O(V + E)` -- optimal for unweighted graphs
- **Space complexity** -- `O(V)` for the queue and distance array
- **Correctness** -- vertices are discovered in non-decreasing distance order; once dequeued, distance is final

> BFS is the baseline. When edges have varying weights, we need more sophisticated approaches -- but for unit-weight graphs, nothing beats BFS.

---

## Slide 04 -- Dijkstra's Algorithm

Dijkstra's algorithm (1959) solves single-source shortest paths for graphs with **non-negative edge weights**.

### Core idea

Greedily extract the vertex with the smallest tentative distance. Its distance is then final (proven optimal). Relax all its outgoing edges.

```
Dijkstra(G, w, s):
  for each vertex v in V:
    dist[v] = inf
    prev[v] = NIL
  dist[s] = 0
  Q = priority queue with all vertices keyed by dist
  while Q is not empty:
    u = extract-min(Q)
    for each edge (u, v) in E:
      if dist[u] + w(u,v) < dist[v]:
        dist[v] = dist[u] + w(u,v)
        prev[v] = u
        decrease-key(Q, v, dist[v])
```

### Key constraint

All edge weights must be non-negative: `w(u, v) >= 0` for all `(u, v) in E`. Negative weights break the greedy guarantee.

---

## Slide 05 -- Dijkstra -- Priority Queue Implementation

The choice of priority queue determines Dijkstra's time complexity:

| Priority Queue | extract-min | decrease-key | Total Complexity |
|----------------|------------|-------------|-----------------|
| **Array (unsorted)** | `O(V)` | `O(1)` | `O(V^2)` |
| **Binary heap** | `O(log V)` | `O(log V)` | `O((V + E) log V)` |
| **Fibonacci heap** | `O(log V)` amortised | `O(1)` amortised | `O(V log V + E)` |

### Practical notes

- **Binary heap** is the standard choice -- simple, cache-friendly, fast in practice
- **Fibonacci heap** wins asymptotically for dense graphs (`E >> V`) but has high constant factors
- **Lazy deletion** variant: push new `(dist, v)` pairs without decrease-key; skip stale entries on extraction. Simpler to implement, same asymptotic cost with binary heap.

```python
import heapq

def dijkstra(graph, source):
    dist = {v: float('inf') for v in graph}
    dist[source] = 0
    pq = [(0, source)]
    while pq:
        d, u = heapq.heappop(pq)
        if d > dist[u]:
            continue          # stale entry
        for v, w in graph[u]:
            if dist[u] + w < dist[v]:
                dist[v] = dist[u] + w
                heapq.heappush(pq, (dist[v], v))
    return dist
```

---

## Slide 06 -- Dijkstra -- Proof of Correctness

### Theorem

When Dijkstra extracts vertex `u` from the priority queue, `dist[u] = delta(s, u)`.

### Proof sketch (by contradiction)

1. Suppose `u` is the first vertex extracted with `dist[u] > delta(s, u)`
2. Let `p` be a true shortest path `s ~> u`. Let `(x, y)` be the first edge on `p` that crosses from the settled set `S` to the unsettled set `V - S`
3. Since `x in S`, `dist[x] = delta(s, x)` (by assumption, `u` is the *first* violation)
4. After relaxing `(x, y)`: `dist[y] <= dist[x] + w(x, y) = delta(s, y)`
5. Since all weights are non-negative: `delta(s, y) <= delta(s, u) < dist[u]`
6. But then `dist[y] < dist[u]`, so the priority queue would have extracted `y` before `u` -- contradiction

> The non-negative weight assumption is essential in step 5. With negative edges, `delta(s, y)` could exceed `delta(s, u)` via a later negative edge, and the proof breaks.

---

## Slide 07 -- Bellman-Ford Algorithm

Bellman-Ford handles graphs with **negative edge weights** (but no negative cycles reachable from `s`).

```
Bellman-Ford(G, w, s):
  for each vertex v in V:
    dist[v] = inf
  dist[s] = 0

  for i = 1 to |V| - 1:       // V-1 passes
    for each edge (u, v) in E:
      if dist[u] + w(u,v) < dist[v]:
        dist[v] = dist[u] + w(u,v)
        prev[v] = u

  // Negative cycle check
  for each edge (u, v) in E:
    if dist[u] + w(u,v) < dist[v]:
      return "Negative cycle detected"
```

### Why `|V| - 1` passes?

A shortest path in a graph with `|V|` vertices has at most `|V| - 1` edges. After pass `i`, the algorithm has found all shortest paths using at most `i` edges.

| Property | Value |
|----------|-------|
| **Time complexity** | `O(V * E)` |
| **Space complexity** | `O(V)` |
| **Handles negative weights** | Yes |
| **Detects negative cycles** | Yes |

---

## Slide 08 -- Bellman-Ford -- Negative Cycle Detection

### What is a negative cycle?

A cycle whose total edge weight is negative. If reachable from `s`, shortest-path distances are `-inf` for all vertices reachable from the cycle -- you can loop forever reducing cost.

```
Example negative cycle:

  A --( 2)--> B --(-5)--> C --( 1)--> A
  
  Cycle weight: 2 + (-5) + 1 = -2
  Each traversal reduces total path weight by 2
```

### Detection mechanism

After `|V| - 1` relaxation passes, all shortest paths are finalized (if no negative cycle). A `|V|`-th pass that still finds improvements proves a negative cycle exists.

### Finding the cycle

- Run the extra pass; record any vertex `v` that improves
- Follow `prev` pointers from `v` for `|V|` steps to ensure you are on the cycle
- Trace `prev` pointers until you revisit a vertex -- that is the cycle

### Applications

- **Currency arbitrage** -- model exchange rates as `-log(rate)` edge weights; a negative cycle means profit
- **Network routing** -- detect routing loops in distance-vector protocols

---

## Slide 09 -- SPFA -- Shortest Path Faster Algorithm

SPFA is a queue-based optimisation of Bellman-Ford. Instead of relaxing all edges in every pass, it only re-relaxes vertices whose distances have changed.

```
SPFA(G, w, s):
  for each vertex v in V:
    dist[v] = inf
    in_queue[v] = false
  dist[s] = 0
  Q = empty queue
  enqueue(Q, s)
  in_queue[s] = true

  while Q is not empty:
    u = dequeue(Q)
    in_queue[u] = false
    for each edge (u, v) in E:
      if dist[u] + w(u,v) < dist[v]:
        dist[v] = dist[u] + w(u,v)
        prev[v] = u
        if not in_queue[v]:
          enqueue(Q, v)
          in_queue[v] = true
```

### Properties

| Property | Value |
|----------|-------|
| **Average-case** | `O(E)` on random graphs |
| **Worst-case** | `O(V * E)` -- same as Bellman-Ford |
| **Handles negative weights** | Yes |
| **Negative cycle detection** | Count relaxations per vertex; if any vertex is relaxed `>= V` times, a negative cycle exists |

> SPFA is popular in competitive programming for its simplicity and good average-case performance. In production systems, Dijkstra with a binary heap is preferred for non-negative weights.

---

## Slide 10 -- DAG Shortest Paths via Topological Sort

For directed acyclic graphs (DAGs), we can find shortest paths in `O(V + E)` by processing vertices in topological order -- even with negative edge weights.

```
DAG-Shortest-Paths(G, w, s):
  topological-sort(G)
  for each vertex v in V:
    dist[v] = inf
  dist[s] = 0

  for each vertex u in topological order:
    for each edge (u, v) in E:
      if dist[u] + w(u,v) < dist[v]:
        dist[v] = dist[u] + w(u,v)
        prev[v] = u
```

### Why it works

In topological order, all predecessors of `v` are processed before `v`. When we relax edges from `u`, we know `dist[u]` is already optimal.

### Properties

- **Time** -- `O(V + E)` -- linear, optimal
- **Handles negative weights** -- yes, because there are no cycles at all
- **No priority queue needed** -- topological order replaces it

### Applications

- **Critical path analysis** -- negate weights, find shortest path = longest path in original
- **Project scheduling** -- PERT/CPM networks are DAGs
- **Dynamic programming** -- many DP problems are shortest paths in an implicit DAG

---

## Slide 11 -- Floyd-Warshall -- All-Pairs Shortest Paths

Floyd-Warshall computes shortest paths between **every pair** of vertices using dynamic programming over intermediate vertices.

```
Floyd-Warshall(G, w):
  let dist[1..V][1..V] = inf
  for each edge (u, v) with weight w:
    dist[u][v] = w
  for each vertex v:
    dist[v][v] = 0

  for k = 1 to V:
    for i = 1 to V:
      for j = 1 to V:
        if dist[i][k] + dist[k][j] < dist[i][j]:
          dist[i][j] = dist[i][k] + dist[k][j]
```

### Recurrence

`dist_k(i, j) = min( dist_{k-1}(i, j), dist_{k-1}(i, k) + dist_{k-1}(k, j) )`

At step `k`, we consider whether routing through vertex `k` improves the path from `i` to `j`.

| Property | Value |
|----------|-------|
| **Time** | `O(V^3)` |
| **Space** | `O(V^2)` |
| **Handles negative weights** | Yes (no negative cycles) |
| **Negative cycle detection** | `dist[v][v] < 0` for some `v` |

> Floyd-Warshall is ideal for dense graphs with small `V`. For sparse graphs, running Dijkstra from every vertex is faster: `O(V * E log V)`.

---

## Slide 12 -- Johnson's Algorithm

Johnson's algorithm computes all-pairs shortest paths efficiently on **sparse graphs** by combining Bellman-Ford and Dijkstra.

### Strategy

1. Add a new vertex `q` with zero-weight edges to all other vertices
2. Run Bellman-Ford from `q` to compute `h(v)` -- a potential function
3. Reweight all edges: `w'(u, v) = w(u, v) + h(u) - h(v)` -- guaranteed non-negative
4. Remove `q` and run Dijkstra from every vertex using reweighted edges
5. Recover original distances: `delta(u, v) = delta'(u, v) - h(u) + h(v)`

### Why reweighting works

For any path `p` from `u` to `v`, the reweighted path weight is `w'(p) = w(p) + h(u) - h(v)`. The `h` values telescope, so the relative ordering of paths is preserved. Shortest paths under `w` remain shortest under `w'`.

| Property | Value |
|----------|-------|
| **Time** | `O(V * E log V)` with binary heap Dijkstra |
| **Space** | `O(V^2)` for the distance matrix |
| **Better than Floyd-Warshall when** | `E = o(V^2)` -- sparse graphs |

> Johnson's is the standard choice for all-pairs on sparse graphs. Floyd-Warshall wins on dense graphs due to simpler constant factors and cache behaviour.

---

## Slide 13 -- A* Search

A* is a best-first search for single-pair shortest paths that uses a **heuristic** to focus exploration toward the goal, dramatically reducing the search space.

### Core idea

Prioritise vertices by `f(v) = g(v) + h(v)`:

- `g(v)` -- actual cost from source to `v` (same as Dijkstra's `dist[v]`)
- `h(v)` -- estimated cost from `v` to goal (the heuristic)
- `f(v)` -- estimated total cost of the cheapest path through `v`

```
A-Star(G, w, s, goal):
  open = priority queue keyed by f
  g[s] = 0
  f[s] = h(s)
  insert(open, s)

  while open is not empty:
    u = extract-min(open)
    if u == goal:
      return reconstruct_path(prev, goal)
    for each edge (u, v) in E:
      tentative = g[u] + w(u, v)
      if tentative < g[v]:
        g[v] = tentative
        f[v] = g[v] + h(v)
        prev[v] = u
        insert-or-update(open, v)
```

> When `h(v) = 0` for all `v`, A* degenerates to Dijkstra's algorithm.

---

## Slide 14 -- A* -- Heuristics & Admissibility

### Admissible heuristic

A heuristic `h` is **admissible** if it never overestimates: `h(v) <= delta(v, goal)` for all `v`. Admissibility guarantees A* finds the optimal path.

### Consistent (monotone) heuristic

`h` is **consistent** if `h(u) <= w(u, v) + h(v)` for every edge `(u, v)`. Consistency implies admissibility and ensures each vertex is expanded at most once.

### Common heuristics for grid/map pathfinding

| Heuristic | Formula | When to use |
|-----------|---------|------------|
| **Manhattan** | `|x1-x2| + |y1-y2|` | 4-directional grid |
| **Euclidean** | `sqrt((x1-x2)^2 + (y1-y2)^2)` | Any-angle movement |
| **Chebyshev** | `max(|x1-x2|, |y1-y2|)` | 8-directional grid |
| **Octile** | `max(dx, dy) + (sqrt(2)-1) * min(dx, dy)` | 8-directional, exact |

### The heuristic trade-off

- `h = 0` -- Dijkstra; explores everything; always optimal
- `h = delta` -- perfect heuristic; A* goes straight to the goal; no wasted work
- `h > delta` -- inadmissible; A* is faster but may return sub-optimal paths
- Stronger admissible heuristics explore fewer vertices but cost more to compute

---

## Slide 15 -- Bidirectional Search

Run two simultaneous searches -- one forward from `s`, one backward from `t` -- and terminate when they meet.

### For unweighted graphs (bidirectional BFS)

- Alternate layers between forward and backward BFS
- Terminate when a vertex appears in both visited sets
- Shortest path passes through the meeting vertex

### For weighted graphs (bidirectional Dijkstra)

- Run forward and backward Dijkstra simultaneously
- Alternate extract-min between the two queues
- Maintain best known path `mu = min(dist_f[v] + dist_b[v])` over all `v`
- Terminate when `dist_f[u] + dist_b[u'] >= mu` (both fronts exceed best path)

### Performance

| Metric | Unidirectional | Bidirectional |
|--------|---------------|--------------|
| **BFS explored** | `O(b^d)` | `O(b^{d/2})` |
| **Speedup** | baseline | ~ `b^{d/2}` |

Where `b` = branching factor, `d` = shortest path length. The quadratic reduction in explored vertices makes bidirectional search valuable for large graphs.

> Bidirectional search requires knowing the target vertex in advance and being able to traverse edges in reverse.

---

## Slide 16 -- Contraction Hierarchies

Contraction hierarchies (CH) are a preprocessing-based speed-up technique for shortest-path queries on static graphs -- the backbone of modern road routing engines.

### Preprocessing

1. Order vertices by "importance" (heuristic: edge difference, contracted neighbours, ...)
2. Contract vertices from least to most important
3. When contracting vertex `v`, add **shortcut edges** between its neighbours if the shortest path between them goes through `v`

### Query

- Run bidirectional Dijkstra, but only relax edges going **upward** in the hierarchy
- Forward search goes up from `s`; backward search goes up from `t`
- The meeting point is the highest-importance vertex on the shortest path

### Performance

| Phase | Time |
|-------|------|
| **Preprocessing** | Minutes to hours (one-time, offline) |
| **Query** | Microseconds -- typically explores < 1000 vertices on continental road networks |

> Google Maps, OSRM, and GraphHopper all use contraction hierarchies or related techniques (ALT, transit nodes). CH reduces query time by 3--6 orders of magnitude compared to plain Dijkstra.

---

## Slide 17 -- Applications -- Routing & Network Protocols

### Road network routing

GPS navigation: preprocess the road graph with contraction hierarchies or A* with landmark heuristics. Real-time traffic requires dynamic edge weights and partial re-computation.

### Network routing protocols

| Protocol | Algorithm | Description |
|----------|-----------|------------|
| **OSPF** | Dijkstra | Link-state: each router has the full topology; runs Dijkstra locally |
| **IS-IS** | Dijkstra | Similar to OSPF; used in large ISP backbones |
| **RIP** | Bellman-Ford | Distance-vector: routers share distance tables with neighbours |
| **BGP** | Path-vector | Policy-based; not strictly shortest-path but related |

### Telecommunications

- Minimum-cost routing in telephone switching networks
- Optical network path computation (GMPLS)
- Fibre-optic cable laying -- shortest path through geographic terrain

> The Internet literally runs on shortest-path algorithms. Every packet you send traverses routes computed by Dijkstra (OSPF) or Bellman-Ford (RIP) at the router level.

---

## Slide 18 -- Applications -- Game AI & Beyond

### Game AI pathfinding

- **Grid-based games** -- A* with Manhattan/octile heuristic is the industry standard
- **Navigation meshes** -- A* on a polygon mesh for 3D environments; funnel algorithm for path smoothing
- **Hierarchical pathfinding** -- HPA* (Hierarchical Pathfinding A*): coarse search between regions, fine search within
- **Flow fields** -- for large numbers of units moving to the same destination (RTS games)

### Robotics

- Motion planning: A* or RRT (Rapidly-exploring Random Trees) in configuration space
- Autonomous vehicles: hybrid A* for kinematic constraints

### Social networks

- Shortest path = degrees of separation
- Betweenness centrality (shortest paths through a node) identifies influential users

### Bioinformatics

- Sequence alignment via edit-distance graphs
- Metabolic pathway analysis -- shortest paths in reaction networks

### Operations research

- Vehicle routing, supply chain optimisation, airline crew scheduling

---

## Slide 19 -- Complexity Comparison

| Algorithm | Time | Space | Negative weights | Negative cycle detection | Use case |
|-----------|------|-------|-----------------|------------------------|----------|
| **BFS** | `O(V + E)` | `O(V)` | No | N/A | Unweighted graphs |
| **Dijkstra (binary heap)** | `O((V+E) log V)` | `O(V)` | No | No | Non-negative weights (default) |
| **Dijkstra (Fibonacci)** | `O(V log V + E)` | `O(V)` | No | No | Dense graphs, non-negative |
| **Bellman-Ford** | `O(V * E)` | `O(V)` | Yes | Yes | Negative weights, cycle detection |
| **SPFA** | `O(E)` avg, `O(VE)` worst | `O(V)` | Yes | Yes | Competitive programming |
| **DAG relaxation** | `O(V + E)` | `O(V)` | Yes | N/A (no cycles) | DAGs, critical path |
| **Floyd-Warshall** | `O(V^3)` | `O(V^2)` | Yes | Yes | All-pairs, dense graphs |
| **Johnson's** | `O(V E log V)` | `O(V^2)` | Yes | Yes | All-pairs, sparse graphs |
| **A*** | `O(b^d)` worst | `O(b^d)` | No | No | Single-pair with heuristic |

> **Decision tree:** Unweighted? BFS. DAG? Topological sort. Non-negative weights? Dijkstra. Negative weights? Bellman-Ford. All-pairs? Floyd-Warshall (dense) or Johnson's (sparse). Single-pair with domain knowledge? A*.

---

## Slide 20 -- Summary & Further Reading

### Key takeaways

- Shortest-path problems have optimal sub-structure -- every sub-path of a shortest path is itself shortest
- Dijkstra is the workhorse for non-negative weights; binary heap implementation is the practical default
- Bellman-Ford handles negative weights at the cost of `O(VE)` time; its extra pass detects negative cycles
- Floyd-Warshall solves all-pairs in `O(V^3)` -- ideal for dense graphs with small vertex counts
- Johnson's algorithm beats Floyd-Warshall on sparse graphs by reweighting edges to eliminate negative weights
- A* focuses search with a heuristic; admissibility guarantees optimality
- Contraction hierarchies make continental-scale routing feasible in microseconds
- The choice of algorithm depends on graph structure (DAG, sparse, dense), weight constraints (non-negative, negative), and query type (single-source, all-pairs, single-pair)

### Recommended reading

| Source | Description |
|--------|------------|
| **CLRS** | *Introduction to Algorithms* -- chapters 22-25 cover BFS, Dijkstra, Bellman-Ford, Floyd-Warshall, Johnson's |
| **Kleinberg & Tardos** | *Algorithm Design* -- excellent exposition of shortest-path problems and network flow |
| **Hart, Nilsson, Raphael** | "A Formal Basis for the Heuristic Determination of Minimum Cost Paths" (1968) -- the original A* paper |
| **Geisberger et al.** | "Contraction Hierarchies: Faster and Simpler Hierarchical Routing" (2008) |
| **Sedgewick & Wayne** | *Algorithms* 4th ed. -- clear Java implementations with visualisations |
| **Jeff Erickson** | *Algorithms* (free online) -- [jeffe.cs.illinois.edu/teaching/algorithms](http://jeffe.cs.illinois.edu/teaching/algorithms/) |
