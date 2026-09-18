# Kinetic Influence Audit

A weighted-network analysis and reconstruction project focused on community structure, node influence, inter-community connectivity, shortest-path relationships, and structural resilience.

The project provides a reproducible pipeline for auditing how influence is distributed through a graph and how the network responds when structurally important nodes are removed.

## Overview

Network influence is rarely captured by a single measurement.

A highly connected node is not necessarily a critical bridge, and a node with high betweenness centrality may interact with a neighborhood very differently from another node with a similar score.

Kinetic Influence Audit combines several complementary analyses:

- weighted graph statistics
- betweenness centrality
- community detection
- community-leader identification
- weighted shortest paths
- node influence scoring
- inter-community connectivity
- resilience under targeted node removal

Together, these provide a more complete picture of network structure than a simple degree ranking.

## Repository Structure

```text
kinetic-influence-audit/
├── cheat/
├── environment/
├── solution/
├── tests/
├── instruction.md
└── task.toml
```

### `instruction.md`

Defines the analysis task, expected behavior, supplied evidence, and required output.

### `environment/`

Contains the reproducible execution environment and task data.

### `solution/`

Contains the reference implementation.

### `tests/`

Contains correctness and verification logic.

### `task.toml`

Defines task metadata and execution configuration.

## Analysis Pipeline

Conceptually, the project follows this workflow:

```text
weighted graph
      ↓
basic graph statistics
      ↓
weighted centrality
      ↓
community detection
      ↓
community leaders
      ↓
leader-to-leader paths
      ↓
influence scoring
      ↓
community connectivity
      ↓
targeted-removal resilience
```

Each stage describes a different aspect of the same network.

## Weighted Centrality

Betweenness centrality measures how frequently a node lies on important paths through the graph.

For weighted networks, shortest paths must account for edge weights rather than treating every edge as equivalent.

This makes centrality useful for identifying structural bridges and nodes that connect otherwise distant regions of the graph.

## Community Detection

The graph is partitioned into communities using a modularity-oriented community-detection method.

Each node is associated with a community, allowing later analysis to distinguish:

- local connectivity
- cross-community connectivity
- community leaders
- inter-community paths

Community structure provides a higher-level view of the graph beyond individual nodes and edges.

## Community Leaders

For each detected community, the pipeline identifies a representative leader using betweenness centrality.

Conceptually:

```text
community
    ↓
member centralities
    ↓
highest-centrality member
    ↓
community leader
```

These leaders provide anchor nodes for examining communication and routing between communities.

## Weighted Leader Paths

Shortest paths are calculated between community leaders using weighted graph distances.

The analysis records information such as:

- path endpoints
- number of traversed edges
- weighted path length
- full node path

This reveals how major regions of the network connect to one another.

## Influence Score

Node influence combines global and local network structure.

The core idea is:

```text
influence
    =
betweenness centrality
    × normalized node degree
    × normalized average neighbor degree
```

This means a node receives a strong influence score when it is simultaneously:

- structurally important to paths through the network
- well connected itself
- surrounded by well-connected neighbors

The measure therefore captures more than raw degree or centrality alone.

## Inter-Community Connectivity

The project audits every pair of communities and measures the links crossing between them.

Relevant metrics include:

```text
cross-community edge count
average cross-community edge weight
```

This helps identify communities that are:

- strongly coupled
- weakly connected
- dependent on a small number of bridges
- relatively isolated

## Network Resilience

Influential nodes can also represent structural vulnerabilities.

The resilience analysis progressively removes high-centrality nodes and measures how the remaining network fragments.

A typical process is:

```text
original graph
     ↓
remove highest-centrality node
     ↓
measure largest connected component
     ↓
remove next node
     ↓
measure again
     ↓
continue
```

This produces a degradation profile showing how quickly connectivity deteriorates under targeted failures.

## Why Multiple Metrics Matter

Consider two nodes with the same degree.

One might belong to a dense local cluster, while the other might be the only bridge connecting two major communities.

Degree alone would treat them similarly.

Combining centrality, neighborhood structure, community information, path analysis, and resilience reveals their different structural roles.

## Technical Areas

This project exercises:

- graph theory
- weighted graphs
- network analysis
- betweenness centrality
- Louvain community detection
- Dijkstra shortest paths
- influence scoring
- graph resilience
- connected components
- algorithmic optimization
- deterministic data processing
- Python scientific computing

## Correctness

Optimized implementations should preserve the semantics of the original analysis.

In particular, performance improvements should not change:

- weighted shortest-path behavior
- community assignments expected by the task
- leader selection
- influence calculations
- cross-community statistics
- resilience measurements

The objective is to improve or reconstruct the analysis while retaining equivalent observable results.

## Goal

Kinetic Influence Audit provides a structured way to examine not only which nodes appear influential, but why they matter to the larger network.

By combining community structure, weighted routing, local neighborhood properties, and targeted-failure analysis, the project produces a broader audit of influence and structural robustness.

## License

No license is currently specified.
