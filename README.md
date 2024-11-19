# Computing approximations of projections

This repository contains an implementation of a method ("sample-and-separate") for computing approximations of projections of polytopes presented in the paper

> Capacity Representation in Sales and Operations Planning: Aggregation through Projection (2024)
> 

The main functionality is implemented in the file `projection.jl`.

## Usage

The implementation is based on the Julia [JuMP](https://jump.dev/) package and requires a linear programming solver to be installed. Within this example, we use the [Gurobi](https://www.gurobi.com/) solver.

```julia
using Gurobi
solver = Gurobi
```

Suppose that the polytope P of interest is given via an extended formulation by $P = \{Tx + t : Ax \le b\}$, where $T$ and $A$ are matrices, $t$ and $b$ are column vectors, of appropriate dimensions. Below, we depict these objects according to Example 1 and Example 2 within the paper.

```julia
# Example 1
A = [1 0 4 0 1 0 0; 0 0 0 5 0 2 0; 0 6 0 0 0 0 3; -1 0 0 0 0 0 0; 0 -1 0 0 0 0 0; 0 0 -1 0 0 0 0; 0 0 0 -1 0 0 0; 0 0 0 0 -1 0 0; 0 0 0 0 0 -1 0; 0 0 0 0 0 0 -1]
b = [10, 20, 30, 0, 0, 0, 0, 0, 0, 0]
T = [1 1 0 0 0 0 0; 0 0 1 1 0 0 0; 0 0 0 0 1 1 1]
t = [0, 0, 0]
```

```julia
# Example 2
A = [6.0 6.0 10.0 0.0; 0.0 0.0 0.0 4.0; 15.0 0.0 7.0 0.0; 0.0 10.0 0.0 9.0; -1.0 0.0 0.0 0.0; 0.0 -1.0 0.0 0.0; 0.0 0.0 -1.0 0.0; 0.0 0.0 0.0 -1.0]
b = [100.0, 100.0, 100.0, 100.0, 0.0, 0.0, 0.0, 0.0]
T = [1.0 1.0 0.0 0.0; 0.0 0.0 1.0 1.0]
t = [0.0, 0.0]
```

In order to replicate one of the instances from experiments in Section 4 of the paper, please consult the data contained in the `experiments/*` folders. As an example, instance 3 of table 2 can be included via

```julia
include("experiments/single-stage/table-2/instance-03.jl")
```

As explained in the paper, the following two parameters should be specified in the algorithm.

```julia
threshold = 0.01 # called θ in the paper
termcount = 100 # called K in the paper
```

An approximation can be computed as follows:

```julia
approximation = Projection.approximate_projection(A, b, T, t, threshold, termcount, solver)
display(approximation)
```
