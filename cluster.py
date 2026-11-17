"""Multiplet building and relative alignment of repeating earthquakes.

Part of inner-core-repeaters by MV2, Munim, Inc. (made with AI agents). Apache-2.0.

Edge convention: edges[(i, j)] = (cc, lag) means event j's arrivals are delayed by
`lag` samples relative to event i's (in their predicted-time frames).
"""


class UnionFind:
    def __init__(self):
        self.p = {}

    def find(self, a):
        self.p.setdefault(a, a)
        root = a
        while self.p[root] != root:
            root = self.p[root]
        while self.p[a] != root:
            self.p[a], a = root, self.p[a]
        return root

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[rb] = ra
        return ra != rb


def components(nodes, edges):
    """Connected components (lists, each sorted) of the graph given by edge keys."""
    uf = UnionFind()
    for n in nodes:
        uf.find(n)
    for (i, j) in edges:
        uf.union(i, j)
    groups = {}
    for n in nodes:
        groups.setdefault(uf.find(n), []).append(n)
    return [sorted(g) for g in groups.values()]


def max_spanning_offsets(members, edges, root):
    """Offsets (samples) of each member relative to root along a maximum-cc spanning tree."""
    mem = set(members)
    es = sorted(((cc, i, j, lag) for (i, j), (cc, lag) in edges.items() if i in mem and j in mem),
                key=lambda e: -e[0])
    uf = UnionFind()
    adj = {m: [] for m in members}
    for cc, i, j, lag in es:
        if uf.union(i, j):
            adj[i].append((j, lag))
            adj[j].append((i, -lag))
    off = {root: 0.0}
    stack = [root]
    while stack:
        u = stack.pop()
        for v, l in adj[u]:
            if v not in off:
                off[v] = off[u] + l
                stack.append(v)
    return off
