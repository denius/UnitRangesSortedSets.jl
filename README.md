# UnitRangesSortedSets

[![Build Status](https://github.com/denius/UnitRangesSortedSets.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/denius/UnitRangesSortedSets.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/denius/UnitRangesSortedSets.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/denius/UnitRangesSortedSets.jl)


Sorted set of `UnitRange`s. Sorted in ascending order and no one range overlaps with another.

    mutable struct UnitRangesSortedSet{K, TU} <: AbstractSet{TU}

`UnitRangesSortedSet` can be created like the standard `Set`:

    urs = UnitRangesSortedSet(somecontainer)

or with `push!`:

```julia-repl
julia> urs = UnitRangesSortedSet{Int}()
UnitRangesSortedSet{Int64}()

julia> push!(urs, 1)
UnitRangesSortedSet{Int64}():
  1:1

julia> push!(urs, 2)
UnitRangesSortedSet{Int64}():
  1:2

julia> push!(urs, 10:12)
UnitRangesSortedSet{Int64}():
   1:2
  10:12
```

Iterating

```julia-repl
julia> for r in urs @show(r) end
r = 1:2
r = 10:12

julia> for r in urs, i in r @show(i) end
i = 1
i = 2
i = 10
i = 11
i = 12

julia> for i in Iterators.flatten(urs) @show(i) end
i = 1
i = 2
i = 10
i = 11
i = 12
```

Deleting elements and ranges:
```julia-repl
julia> delete!(urs, 10:11)
UnitRangesSortedSet{Int64}():
   1:2
  12:12

julia> delete!(urs, 1)
UnitRangesSortedSet{Int64}():
   2:2
  12:12
```


Some internal details:
The struct contains `ranges::SortedDict{K,K}`, where the key is `first(range)`, and the value is `last(range)`.
Note: for `Char` the `StepRange{Char,UInt8}` with `oneunit(UInt8)` step will be created.
