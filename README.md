# UnitRangesSortedSets

[![Build Status](https://github.com/denius/UnitRangesSortedSets.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/denius/UnitRangesSortedSets.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/denius/UnitRangesSortedSets.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/denius/UnitRangesSortedSets.jl)


Sorted set of `UnitRange`s. Sorted in ascending order and no one range overlaps with another.

    mutable struct UnitRangesSortedSet{K, TU} <: AbstractSet{TU}

`UnitRangesSortedSet` can be created like the standard `Set`:

    UnitRangesSortedSet(somecontainer)

```julia
julia> using UnitRangesSortedSets

julia> UnitRangesSortedSet((1, 2, 4))
UnitRangesSortedSet{Int64}():
  1:2
  4:4

julia> UnitRangesSortedSet(('a':'z', 'α':'ω'))
UnitRangesSortedSet{Char}():
  'a':'z'
  'α':'ω'

julia> Random.seed!(1234);

julia> UnitRangesSortedSet(rand(1:20, 10))
UnitRangesSortedSet{Int64}():
   5:5
   7:8
  10:11
  15:16
  18:18
  20:20
```

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

Iterating over set of ranges:

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

julia> collect(urs)
2-element Vector{UnitRange{Int64}}:
 1:2
 10:12
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

# SubSet

It is possible to create subsets of `UnitRangesSortedSet` (like a `view` for `Array`s):
```julia
julia> urs = UnitRangesSortedSet((1:2, 10:12))
UnitRangesSortedSet{Int64}():
   1:2
  10:12

julia> ss = subset(urs, 0:10)
2-element subset(UnitRangesSortedSet{Int64}, DataStructures.Tokens.IntSemiToken(3):DataStructures.Tokens.IntSemiToken(4)):
   1:2
  10:10
```

The `subset` object is a static, iterable view of the container.

# Two types of `UnitRangesSortedSet`

The first type `UnitRangesSortedSet{K}` contains `SortedDict{K,K}`,
```julia
mutable struct UnitRangesSortedSet{K,TU} <: AbstractUnitRangesSortedContainer{K,TU}
    ranges::SortedDict{K,K,FOrd}
end
```
where each element of the dict contains the `first(range)` as key, and the `last(range)` as value.

The second implementation `PlainUnitRangesSortedSet{K}` is based on `Vector{K}`s:
```julia
mutable struct PlainUnitRangesSortedSet{K,TU} <: AbstractUnitRangesSortedContainer{K,TU}
    rstarts::Vector{K}
    rstops::Vector{K}
end
```
where `rstarts::Vector{K}` and `rstops::Vector{K}` are the starts and stops of
the ranges respectively.

These two implementations have a similar API but different speeds.
When created from elements in random order, `UnitRangesSortedSet` is vastly superior
to the `Plain` variant.
But in searching operations (`in()`, `subset()`) `PlainUnitRangesSortedSet` variant is faster:
in Julia-v1.6 it is twice as fast, in Julia-1.8 the speedup is about 20-30%.
In either case, both of them can be converted to each other using the appropriate constructor.

Note: for `Char`, `StepRange{Char,UInt8}` will be used, with a step of `oneunit(UInt8)` if needed.

