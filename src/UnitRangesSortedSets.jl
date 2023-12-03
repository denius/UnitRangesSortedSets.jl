
module UnitRangesSortedSets

export AbstractUnitRangesSortedSet, AbstractUnitRangesSortedContainer, AbstractUnitRangesSortedSubSet
export UnitRangesSortedSet, VecUnitRangesSortedSet
export UnitRangesSortedSubSet, UnitRangesSortedSubSet1, UnitRangesSortedSubSet0
export searchsortedrange, searchsortedrangefirst, searchsortedrangelast
export getrange, getindex, beforefirstindex, pastlastindex
export subset


using Base: ForwardOrdering, Forward
const FOrd = ForwardOrdering

using DataStructures
using DataStructures: DataStructures.Tokens.IntSemiToken, DataStructures.SDMToken


safe_add(x::T, y) where {T} = x + T(y)
safe_add(x::T, y) where {T<:Integer} = ((z, flag) = Base.add_with_overflow(x, T(y));   !flag ? z : typemax(T))
safe_add(x::T, y) where {T<:AbstractChar} = (z = UInt64(x) + UInt64(y);   z < UInt128(2)^(8 * sizeof(T)) ? T(z) : typemax(T))

safe_sub(x::T, y) where {T} = x - T(y)
safe_sub(x::T, y) where {T<:Integer} = ((z, flag) = Base.sub_with_overflow(x, T(y));   !flag ? z : typemin(T))
safe_sub(x::T, y) where {T<:AbstractChar} = (z = Int64(x) - Int64(y);   z >= 0 ? T(z) : typemin(T))


@inline to_urange(::Type{TU}, l::L, r::L) where {TU<:UnitRange,L} =
    TU(l, r)
@inline to_urange(::Type{TU}, l::L, r::R) where {TU<:UnitRange,L,R} =
    TU(promote_type(L, R)(l), promote_type(L, R)(r))
@inline to_urange(::Type{TU}, kk::AbstractRange) where {TU<:UnitRange} =
    TU(first(kk), last(kk))
@inline to_urange(::Type{TU}, l, r) where {TU<:StepRange{T,S}} where {T,S} =
    TU(l, oneunit(S), r)
@inline to_urange(::Type{TU}, kk::AbstractRange) where {TU<:StepRange{T,S}} where {T,S} =
    TU(first(kk), oneunit(S), last(kk))


@inline inferrangetype(::Type{K}) where {K<:Real} = UnitRange{K}
@inline function inferrangetype(::Type{K}) where K
    T = typeof(K(0):K(0))
    return T <: StepRange ? StepRange{K,UInt8} : T
end


abstract type AbstractUnitRangesSortedSet{K,TU} <: AbstractSet{TU} end

abstract type AbstractUnitRangesSortedContainer{K,TU} <: AbstractUnitRangesSortedSet{K,TU} end
abstract type AbstractUnitRangesSortedSubSet{K,TU,P} <: AbstractUnitRangesSortedSet{K,TU} end


"""
Sorted set of `UnitRange`s. Sorted in ascending order and no one range overlaps with another.

```julia
mutable struct UnitRangesSortedSet{K,TU} <: AbstractSet{TU}
    "Index of last used range."
    lastusedrangeindex::IntSemiToken
    "Storage for ranges: the key of Dict is the `first(range)`, and the value of Dict is the `last(range)`."
    ranges::SortedDict{K,K,FOrd}
end
```

Ranges stored in `SortedDict` as `first`s in keys and `last`s in values.

`UnitRangesSortedSet` can be created like the standard `Set`:

```julia
    UnitRangesSortedSet(somecontainer)
```

for example
```julia
julia> using UnitRangesSortedSets

julia> UnitRangesSortedSet((1, 2, 4))
UnitRangesSortedSet{Int64} with 2 elements:
  1:2
  4:4

julia> UnitRangesSortedSet(('a':'z', 'α':'ω'))
UnitRangesSortedSet{Char} with 2 elements:
  'a':'z'
  'α':'ω'

julia> Random.seed!(1234);

julia> UnitRangesSortedSet(rand(1:20, 10))
UnitRangesSortedSet{Int64} with 6 elements:
   5:5
   7:8
  10:11
  15:16
  18:18
  20:20
```

or with `push!`:

```julia
julia> urs = UnitRangesSortedSet{Int}()
UnitRangesSortedSet{Int64}()

julia> push!(urs, 1)
UnitRangesSortedSet{Int64} with 1 element:
  1:1

julia> push!(urs, 2)
UnitRangesSortedSet{Int64} with 1 element:
  1:2

julia> push!(urs, 10:12)
UnitRangesSortedSet{Int64} with 2 elements:
   1:2
  10:12
```

Iterating over set of ranges:

```julia
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
```julia
julia> delete!(urs, 10:11)
UnitRangesSortedSet{Int64} with 2 elements:
   1:2
  12:12

julia> delete!(urs, 1)
UnitRangesSortedSet{Int64} with 2 elements:
   2:2
  12:12
```

Note: for `Char` the `StepRange{Char,UInt8}` with `oneunit(UInt8)` step will be created.
"""
mutable struct UnitRangesSortedSet{K,TU} <: AbstractUnitRangesSortedContainer{K,TU}
    "Index of last used range."
    lastusedrangeindex::IntSemiToken
    "Storage for ranges: the key of Dict is the `first(range)`, and the value of Dict is the `last(range)`."
    ranges::SortedDict{K,K,FOrd}
end

function UnitRangesSortedSet{K,TU}() where {K,TU}
    ranges = SortedDict{K,K,FOrd}(Forward)
    UnitRangesSortedSet{K,TU}(beforestartsemitoken(ranges), ranges)
end
UnitRangesSortedSet{K}() where {K} = UnitRangesSortedSet{K,inferrangetype(K)}()
function UnitRangesSortedSet(rs::AbstractUnitRangesSortedSet{K,TU}) where {K,TU}
    ranges = SortedDict{K,K,FOrd}(Forward)
    for r in rs
        ranges[first(r)] = last(r)
    end
    UnitRangesSortedSet{K,TU}(beforestartsemitoken(ranges), ranges)
end


"""
Sorted set of `UnitRange`s. Sorted in ascending order and no one range overlaps with another.

```julia
mutable struct VecUnitRangesSortedSet{K,TU} <: AbstractSet{TU}
    "Index of last used range."
    lastusedrangeindex::Int
    "Storage for ranges starts,"
    rstarts::Vector{K}
    "and stops."
    rstops::Vector{K}
end
```

Ranges stored in `Vector`s: in `rstarts` stored `first`s of ranges, and in `rstops` stored `last`s.

For further details see `UnitRangesSortedSet`.
"""
mutable struct VecUnitRangesSortedSet{K,TU} <: AbstractUnitRangesSortedContainer{K,TU}
    "Index of last used range."
    lastusedrangeindex::Int
    "Storage for ranges starts,"
    rstarts::Vector{K}
    "and stops."
    rstops::Vector{K}
end

VecUnitRangesSortedSet{K,TU}() where {K,TU} = VecUnitRangesSortedSet{K,TU}(0, Vector{K}(undef, 0), Vector{K}(undef, 0))
VecUnitRangesSortedSet{K}() where {K} = VecUnitRangesSortedSet{K,inferrangetype(K)}()
function VecUnitRangesSortedSet(rs::AbstractUnitRangesSortedSet{K,TU}) where {K,TU}
    rstarts = Vector{K}(undef, length(rs))
    rstops = Vector{K}(undef, length(rs))
    for (i, r) in enumerate(rs)
        rstarts[i] = first(r)
        rstops[i] = last(r)
    end
    VecUnitRangesSortedSet{K,TU}(firstindex(rstarts) - 1, rstarts, rstops)
end


#
# SubSets
#


"""
It is possible to create the subset of any `AbstractUnitRangesSortedSet`, like a `view` for `Array`s:
```julia
julia> urs = UnitRangesSortedSet((1:2, 10:12))
UnitRangesSortedSet{Int64} with 2 elements:
   1:2
  10:12

julia> ss = subset(urs, 0:10)
2-element subset(UnitRangesSortedSet{Int64}, DataStructures.Tokens.IntSemiToken(3):DataStructures.Tokens.IntSemiToken(4)):
   1:2
  10:10
```

The `subset` object is an static, iterable view of the container

    struct UnitRangesSortedSubSet{K, TU, P, Tix} <: AbstractUnitRangesSortedSubSet{K, TU, P}

"""
struct UnitRangesSortedSubSet{K,TU,P,Tix} <: AbstractUnitRangesSortedSubSet{K,TU,P}
    "The `<:AbstractUnitRangesSortedSet` which subset point to."
    parent::P
    "`kstart:kstop` is the bounds of subset."
    kstart::K
    kstop::K
    "Saved first range `parent[firstindex]` truncated to `kstart` from left."
    firstrange::TU
    "Saved last range `parent[lastindex]` truncated to `kstop` from right."
    lastrange::TU
    "Range `firstindex:lastindex` locates the indices diapason in `parent` where `kstart:kstop` should be inserted."
    firstindex::Tix
    lastindex::Tix
    "Pre `firstindex` position, used for acceleration of `iterate`."
    beforefirstindex::Tix
    "Post `lastindex` position, used for acceleration of `iterate`."
    pastlastindex::Tix
    "Length of subset."
    numranges::Int
end


"""
    struct UnitRangesSortedSubSet1{K, TU, P, Tix} <: AbstractUnitRangesSortedSubSet{K, TU, P}

is like the `UnitRangesSortedSubSet` but with only one range or part of range in subset.
"""
struct UnitRangesSortedSubSet1{K,TU,P,Tix} <: AbstractUnitRangesSortedSubSet{K,TU,P}
    "The `<:AbstractUnitRangesSortedSet` which subset point to."
    parent::P
    "`kstart:kstop` is the bounds of subset."
    kstart::K
    kstop::K
    "Single continuous range `kstart:kstop` from `parent` stored for simplicity."
    singlerange::TU
    "Range `firstindex:lastindex` locates the point in `parent` where `kstart:kstop` should be inserted."
    firstindex::Tix
    lastindex::Tix
    "Pre `firstindex` position, used for acceleration of `iteration`."
    beforefirstindex::Tix
    "Post `lastindex` position, used for acceleration of `iteration`."
    pastlastindex::Tix
    "Length of subset."
    numranges::Int
end


"""
    struct UnitRangesSortedSubSet0{K, TU, P, Tix} <: AbstractUnitRangesSortedSubSet{K, TU, P}

is like the `UnitRangesSortedSubSet` but with the zero-length subset range in subset.
"""
struct UnitRangesSortedSubSet0{K,TU,P,Tix} <: AbstractUnitRangesSortedSubSet{K,TU,P}
    #"Empty range `kstart:kstop` for simplicity."
    #singlerange::TU
    "The `<:AbstractUnitRangesSortedSet` which subset point to."
    parent::P
    "`kstart:kstop` is the bounds of subset."
    kstart::K
    kstop::K
    "Empty range `firstindex:lastindex` locates the point in `parent` where `kstart:kstop` should be inserted."
    firstindex::Tix
    lastindex::Tix
    "Pre `firstindex` position, used for acceleration of `iteration`."
    beforefirstindex::Tix
    "Post `lastindex` position, used for acceleration of `iteration`."
    pastlastindex::Tix
    "Length of subset."
    numranges::Int
end


subset(rs::P, ::Colon) where {P<:AbstractUnitRangesSortedSet{K,TU}} where {K,TU} =
    subset(rs, to_urange(TU, first(first(rs)), last(last(rs))))

subset(rs::P, kk::AbstractRange) where {P<:AbstractUnitRangesSortedSubSet{K,TU}} where {K,TU} =
    subset(rs.parent, to_urange(TU, max(first(first(rs)), K(first(kk))), min(last(last(rs)), K(last(kk)))))

function subset(rs::P, kk::AbstractRange) where {P<:AbstractUnitRangesSortedContainer{K,TU}} where {K,TU}
    kk = to_urange(TU, kk)
    iir_kk = searchsortedrange(rs, kk)
    if length(kk) == 0 || length(iir_kk) == 0
        if length(kk) == 0 && length(iir_kk) == 1
            iir_kk = create_indexrange(rs, first(iir_kk), regress(rs, first(iir_kk)))
        end
        if length(kk) != 0
            kk = to_urange(TU, last(kk), first(kk))
        end
        return UnitRangesSortedSubSet0{K,TU,P,typeof(first(iir_kk))}(rs, first(kk), last(kk), first(iir_kk), last(iir_kk),
                                       regress(rs, first(iir_kk)), advance(rs, last(iir_kk)), length(iir_kk))
    elseif length(iir_kk) == 1
        singlerange = getindex(rs, first(iir_kk))
        if first(singlerange) < first(kk)
            singlerange = to_urange(TU, first(kk), last(singlerange))
        end
        if last(kk) < last(singlerange)
            singlerange = to_urange(TU, first(singlerange), last(kk))
        end
        return UnitRangesSortedSubSet1{K,TU,P,typeof(first(iir_kk))}(rs, first(kk), last(kk), singlerange,
                                       first(iir_kk), last(iir_kk), regress(rs, first(iir_kk)), advance(rs, last(iir_kk)), length(iir_kk))
    else
        firstrange = getindex(rs, first(iir_kk))
        if first(firstrange) < first(kk)
            firstrange = to_urange(TU, first(kk), last(firstrange))
        end
        lastrange = getindex(rs, last(iir_kk))
        if last(kk) < last(lastrange)
            lastrange = to_urange(TU, first(lastrange), last(kk))
        end
        return UnitRangesSortedSubSet{K,TU,P,typeof(first(iir_kk))}(rs, first(kk), last(kk), firstrange, lastrange,
                                      first(iir_kk), last(iir_kk), regress(rs, first(iir_kk)), advance(rs, last(iir_kk)), length(iir_kk))
    end
end

@inline Base.copy(rs::T) where {T<:AbstractUnitRangesSortedSubSet{K,TU}} where {K,TU} =
    intersect(rs.parent, to_urange(TU, rs.kstart, rs.kstop))

@inline Base.parent(rs::T) where {T<:AbstractUnitRangesSortedSubSet{K,TU}} where {K,TU} = rs.parent


#
# Convert
#

(::Type{T})(values::Union{AbstractVector, AbstractSet, Tuple}) where {T<:AbstractUnitRangesSortedContainer} =
    eltype(values) <: AbstractRange ? T{eltype(eltype(values))}(values) : T{eltype(values)}(values)

(::Type{T})(values::Union{AbstractVector, AbstractSet, Tuple}) where {T<:AbstractUnitRangesSortedContainer{K}} where {K} =
    T{inferrangetype(K)}(values)

function (::Type{T})(values::Union{AbstractVector, AbstractSet, Tuple}) where {T<:AbstractUnitRangesSortedContainer{K,TU}} where {K,TU}
    rs = T()
    for r in values
        push!(rs, r)
    end
    rs
end

Base.convert(::Type{<:VecUnitRangesSortedSet{K}}, rs::VecUnitRangesSortedSet{K}) where K = rs
Base.convert(::Type{<:UnitRangesSortedSet{K}}, rs::UnitRangesSortedSet{K}) where K = rs
function Base.convert(::Type{T}, rs::AbstractUnitRangesSortedContainer) where {T<:AbstractVector{Tv}} where {Tv<:AbstractRange}
    V = T(undef, length(rs))
    for (i, r) in enumerate(rs)
        V[i] = to_urange(Tv, r)
    end
    V
end
function Base.convert(::Type{T}, rs::AbstractUnitRangesSortedContainer) where {T<:AbstractVector{Tv}} where {Tv}
    V = T(undef, sum(length(r) for r in rs))
    i = 0
    for r in rs, v in r
        V[i+=1] = Tv(v)
    end
    V
end
function Base.convert(::Type{T}, rs::AbstractUnitRangesSortedContainer{K,TU}) where {T<:AbstractVector,K,TU}
    V = T{K}(undef, sum(length(r) for r in rs))
    i = 0
    for r in rs, v in r
        V[i+=1] = K(v)
    end
    V
end
function Base.convert(::Type{T}, rs::AbstractUnitRangesSortedContainer) where {T<:AbstractSet{Tv}} where {Tv<:AbstractRange}
    S = T()
    for r in rs
        push!(S, to_urange(Tv, r))
    end
    S
end
function Base.convert(::Type{T}, rs::AbstractUnitRangesSortedContainer) where {T<:AbstractSet{Tv}} where {Tv}
    S = T()
    for r in rs, v in r
        push!(S, Tv(v))
    end
    S
end
function Base.convert(::Type{T}, rs::AbstractUnitRangesSortedContainer{K}) where {T<:AbstractSet,K}
    S = T{K}()
    for r in rs, v in r
        push!(S, K(v))
    end
    S
end





"Type `UnitRange` for `DataStructures.Tokens.IntSemiToken`."
struct URSSIndexURange{P,Tix} <: AbstractUnitRange{Tix}
    start::Tix
    stop::Tix
    parent::P
    len::Int
end



"""
Returns `UnitRange` and `UnitRange`-alike ranges which are used for searching and indexing ranges in `AbstractUnitRangesSortedSet`.
Returned by `searchsorted()` function to point sought or missing ranges in given set.
The usual `searchsorted()` function returns just `UnitRange`, but `UnitRangesSortedSet` uses
`DataStructures.Tokens.IntSemiToken` to indexing keys/values and it is not possible to create `UnitRange{IntSemiToken}`,
thus was created custom one named `URSSIndexURange`.
"""
@inline create_indexrange(parent::P, l, r) where {P<:VecUnitRangesSortedSet} = UnitRange{Int}(l, r)
@inline create_indexrange(parent::P, l, r) where {P<:UnitRangesSortedSet} = URSSIndexURange(parent, l, r)


@inline function URSSIndexURange(rs::P, l::Tix, r::Tix) where {P<:VecUnitRangesSortedSet, Tix}
    @boundscheck l < beforefirstindex(rs) || l > pastlastindex(rs) && return throw(BoundsError(rs, l))
    @boundscheck r < beforefirstindex(rs) || r > pastlastindex(rs) && return throw(BoundsError(rs, r))
    if l == pastlastindex(rs) ||
       r == beforefirstindex(rs) ||
       l > r
        len = 0
    else
        len = r - l + 1
    end
    URSSIndexURange{P,Tix}(l, r, rs, len)
end

@inline function URSSIndexURange(rs::P, l::Tix, r::Tix) where {P<:UnitRangesSortedSet, Tix<:IntSemiToken}
    @boundscheck status((rs.ranges, l)) == 0 && return throw(KeyError(l))
    @boundscheck status((rs.ranges, r)) == 0 && return throw(KeyError(r))
    if l == pastlastindex(rs) ||
       r == beforefirstindex(rs) ||
       compare(rs.ranges, l, r) == 1
        len = 0
    else
        len = 0
        i = l
        while compare(rs.ranges, i, r) < 1
            len += 1
            i = advance(rs, i)
        end
    end
    URSSIndexURange{P,Tix}(l, r, rs, len)
end

@inline Base.length(ur::URSSIndexURange) = ur.len
@inline Base.step(rs::URSSIndexURange) = 1

@inline Base.:(==)(l::URSSIndexURange{Pl,Tixl}, r::URSSIndexURange{Pr,Tixr}) where
                   {Pl<:UnitRangesSortedSet,Tixl,Pr<:VecUnitRangesSortedSet,Tixr} = false
@inline Base.:(==)(l::URSSIndexURange{Pl,Tixl}, r::URSSIndexURange{Pr,Tixr}) where
                   {Pl<:VecUnitRangesSortedSet,Tixl,Pr<:UnitRangesSortedSet,Tixr} = false
@inline Base.:(==)(l::URSSIndexURange{Pl,Tixl}, r::URSSIndexURange{Pr,Tixr}) where
                   {Pl<:VecUnitRangesSortedSet,Tixl,Pr<:VecUnitRangesSortedSet,Tixr} =
    l.parent === r.parent &&
    l.start == r.start &&
    l.stop == r.stop
@inline Base.:(==)(l::URSSIndexURange{Pl,Tixl}, r::URSSIndexURange{Pr,Tixr}) where
                   {Pl<:UnitRangesSortedSet,Tixl,Pr<:UnitRangesSortedSet,Tixr} =
    l.parent === r.parent &&
    compare(l.parent.ranges, l.start, r.start) == 0 &&
    compare(l.parent.ranges, l.stop, r.stop) == 0

@inline Base.first(ur::URSSIndexURange) = ur.start
@inline Base.last(ur::URSSIndexURange) = ur.stop
@inline Base.firstindex(ur::URSSIndexURange) = ur.start
@inline Base.lastindex(ur::URSSIndexURange) = ur.stop

#@inline Base.collect(ur::URSSIndexURange{P,Tix}) where {P<:VecUnitRangesSortedSet,Tix} =
#    map(s->getindex(ur.parent, s), ur.start:ur.stop)
#@inline Base.collect(ur::URSSIndexURange{P,Tix}) where {P<:UnitRangesSortedSet,Tix} =
#    map(s->getindex(ur.parent, s), onlysemitokens(inclusive(ur.parent.ranges, ur.start, ur.stop)))

@inline DataStructures.advance(ur::URSSIndexURange, st::IntSemiToken) = advance((ur.parent, st))
@inline DataStructures.regress(ur::URSSIndexURange, st::IntSemiToken) = regress((ur.parent, st))
@inline beforefirstindex(ur::URSSIndexURange) = beforestartsemitoken(ur.parent)
@inline pastlastindex(ur::URSSIndexURange) = pastendsemitoken(ur.parent)

Base.length(rs::VecUnitRangesSortedSet) = length(rs.rstarts)
Base.length(rs::UnitRangesSortedSet) = length(rs.ranges)
Base.length(rs::AbstractUnitRangesSortedSubSet) = rs.numranges
Base.haslength(rs::AbstractUnitRangesSortedSet) = true
Base.hasfastin(rs::AbstractUnitRangesSortedSet) = true
Base.isempty(rs::AbstractUnitRangesSortedSet) = length(rs) == 0
Base.size(rs::AbstractUnitRangesSortedSet) = (length(rs),)
#Base.axes(rs::AbstractUnitRangesSortedSet) = (Base.OneTo(rs.n),)
Base.eltype(::AbstractUnitRangesSortedSet{K,TU}) where {K,TU} = TU
#Base.IndexStyle(::AbstractUnitRangesSortedSet) = IndexLinear()


function Base.collect(::Type{ElType}, rs::AbstractUnitRangesSortedSet{K,TU}) where {ElType,K,TU}
    T = inferrangetype(ElType)
    res = Vector{T}(undef, length(rs))
    i = 0
    for r in rs
        res[i+=1] = to_urange(T, ElType(first(r)), ElType(last(r)))
    end
    return res
end
function Base.collect(rs::AbstractUnitRangesSortedSet{K,TU}) where {K,TU}
    res = Vector{TU}(undef, length(rs))
    i = 0
    for r in rs
        res[i+=1] = to_urange(TU, r)
    end
    return res
end

@inline function indexcompare(rs::VecUnitRangesSortedSet, i, j)
    if i < j
        return -1
    elseif i > j
        return 1
    else
        return 0
    end
end
@inline function indexcompare(rs::AbstractUnitRangesSortedSubSet{K,TU,P}, i, j) where {K,TU,P<:VecUnitRangesSortedSet}
    if i < j
        return -1
    elseif i > j
        return 1
    else
        return 0
    end
end
@inline indexcompare(rs::UnitRangesSortedSet, i, j) = compare(rs.ranges, i, j)
@inline indexcompare(rs::AbstractUnitRangesSortedSubSet{K,TU,P}, i, j) where {K,TU,P<:UnitRangesSortedSet} =
    compare(rs.parent.ranges, i, j)


@inline function getrange(rs::AbstractUnitRangesSortedSet{K,TU}, i) where {K,TU}
    if in(i, rs)
        return to_urange(TU, first(i), last(i))
    else
        return nothing
    end
end

@inline getindex_tuple(rs::VecUnitRangesSortedSet, i) = (rs.rstarts[i], rs.rstops[i])
@inline getindex_tuple(rs::UnitRangesSortedSet, i::IntSemiToken) = tuple(deref((rs.ranges, i))...)
@inline getindex_tuple(rs::UnitRangesSortedSubSet0, i) = (r = getindex(rs, i); tuple(first(r), last(r)))
@inline getindex_tuple(rs::UnitRangesSortedSubSet1, i) = (r = getindex(rs, i); tuple(first(r), last(r)))
@inline getindex_tuple(rs::UnitRangesSortedSubSet, i) = (r = getindex(rs, i); tuple(first(r), last(r)))
@inline Base.getindex(rs::VecUnitRangesSortedSet{K,TU}, i) where {K,TU} =
    to_urange(TU, rs.rstarts[i], rs.rstops[i])
@inline Base.getindex(rs::UnitRangesSortedSet{K,TU}, i::IntSemiToken) where {K,TU} =
    to_urange(TU, deref((rs.ranges, i))...)
@inline Base.getindex(rs::UnitRangesSortedSet{K,TU}, t::SDMToken) where {K,TU} =
    to_urange(TU, deref(t)...)
@inline Base.getindex(ur::URSSIndexURange, i::Integer) = _getindex(ur, i) # ambiguity
@inline Base.getindex(ur::URSSIndexURange, i::Colon) = _getindex(ur, i)
@inline Base.getindex(ur::URSSIndexURange, i::StepRange{T}) where {T<:Integer} = _getindex(ur, i)
@inline Base.getindex(ur::URSSIndexURange, i::AbstractUnitRange{T}) where {T<:Integer} = _getindex(ur, i)
@inline function _getindex(ur::URSSIndexURange, i)
    @boundscheck indexcompare(ur.parent, first(ur), i) != 1 &&
                 indexcompare(ur.parent, i, last(ur)) != 1 || throw(BoundsError(ur, i))
    getindex(ur.parent, i)
end
@inline Base.getindex(rs::UnitRangesSortedSubSet0{K,TU}, i) where {K,TU} =
    (throw(BoundsError(rs)); to_urange(TU, rs.kstart, rs.kstop))
@inline function Base.getindex(rs::UnitRangesSortedSubSet1, i)
    @boundscheck i == rs.firstindex || throw(BoundsError(rs, i))
    rs.singlerange
end
@inline function Base.getindex(rs::UnitRangesSortedSubSet{K,TU,P}, i) where {K,TU,P<:VecUnitRangesSortedSet}
    @boundscheck beforefirstindex(rs) < i < pastlastindex(rs) || throw(BoundsError(rs, i))
    if i == rs.firstindex
        return rs.firstrange
    elseif i < rs.lastindex
        return getindex(rs.parent, i)
    else
        return rs.lastrange
    end
end
@inline function Base.getindex(rs::UnitRangesSortedSubSet, i)
    @boundscheck indexcompare(rs.parent, rs.firstindex, i) != 1 &&
                 indexcompare(rs.parent, i, rs.lastindex) != 1 || throw(BoundsError(rs, i))
    if indexcompare(rs, i, rs.firstindex) == 0
        return rs.firstrange
    elseif indexcompare(rs, i, rs.lastindex) == -1
        return getindex(rs.parent, i)
    else
        return rs.lastrange
    end
end

@inline getindex_rangestart(rs::VecUnitRangesSortedSet, i) = rs.rstarts[i]
@inline getindex_rangestop(rs::VecUnitRangesSortedSet, i) = rs.rstops[i]
@inline getindex_rangestart(rs::UnitRangesSortedSet, i) = deref_key((rs.ranges, i))
@inline getindex_rangestop(rs::UnitRangesSortedSet, i) = deref_value((rs.ranges, i))
@inline getindex_rangestart(rs::AbstractUnitRangesSortedSubSet, i) = first(getindex(rs, i))
@inline getindex_rangestop(rs::AbstractUnitRangesSortedSubSet, i) = last(getindex(rs, i))

@inline Base.firstindex(rs::VecUnitRangesSortedSet) = firstindex(rs.rstarts)
@inline Base.firstindex(rs::UnitRangesSortedSet) = startof(rs.ranges)
@inline Base.firstindex(rs::AbstractUnitRangesSortedSubSet) = rs.firstindex
@inline Base.lastindex(rs::VecUnitRangesSortedSet) = lastindex(rs.rstarts)
@inline Base.lastindex(rs::UnitRangesSortedSet) = lastindex(rs.ranges)
@inline Base.lastindex(rs::AbstractUnitRangesSortedSubSet) = rs.lastindex
"""
    first(rs::AbstractUnitRangesSortedSet)

Returns first range from the set of ranges `rs`.
"""
@inline Base.first(rs::VecUnitRangesSortedSet{K,TU}) where {K,TU} =
    to_urange(TU, rs.rstarts[1], rs.rstops[1])
@inline Base.first(rs::UnitRangesSortedSet{K,TU}) where {K,TU} =
    to_urange(TU, deref((rs.ranges, startof(rs.ranges)))...)
@inline Base.first(rs::UnitRangesSortedSubSet0) = getindex(rs, firstindex(rs))
@inline Base.first(rs::UnitRangesSortedSubSet1) = rs.singlerange
@inline Base.first(rs::UnitRangesSortedSubSet) = getindex(rs, firstindex(rs))
"""
    last(rs::AbstractUnitRangesSortedSet)

Returns last range from the set of ranges `rs`.
"""
@inline Base.last(rs::VecUnitRangesSortedSet{K,TU}) where {K,TU} =
    to_urange(TU, rs.rstarts[end], rs.rstops[end])
@inline Base.last(rs::UnitRangesSortedSet{K,TU}) where {K,TU} =
    to_urange(TU, deref((rs.ranges, lastindex(rs.ranges)))...)
@inline Base.last(rs::UnitRangesSortedSubSet0) = getindex(rs, lastindex(rs))
@inline Base.last(rs::UnitRangesSortedSubSet1) = rs.singlerange
@inline Base.last(rs::UnitRangesSortedSubSet) = getindex(rs, lastindex(rs))

@inline beforefirstindex(rs::VecUnitRangesSortedSet) = firstindex(rs.rstarts) - 1
@inline beforefirstindex(rs::UnitRangesSortedSet) = beforestartsemitoken(rs.ranges)
@inline beforefirstindex(rs::AbstractUnitRangesSortedSubSet) = rs.beforefirstindex
@inline pastlastindex(rs::VecUnitRangesSortedSet) = lastindex(rs.rstarts) + 1
@inline pastlastindex(rs::UnitRangesSortedSet) = pastendsemitoken(rs.ranges)
@inline pastlastindex(rs::AbstractUnitRangesSortedSubSet) = rs.pastlastindex

@inline DataStructures.advance(rs::VecUnitRangesSortedSet, state) = state + 1
@inline DataStructures.advance(rs::UnitRangesSortedSet, state) = advance((rs.ranges, state))
@inline DataStructures.advance(rs::AbstractUnitRangesSortedSubSet, state) = advance(rs.parent, state)
@inline DataStructures.regress(rs::VecUnitRangesSortedSet, state) = state - 1
@inline DataStructures.regress(rs::UnitRangesSortedSet, state) = regress((rs.ranges, state))
@inline DataStructures.regress(rs::AbstractUnitRangesSortedSubSet, state) = regress(rs.parent, state)

#
# Searching functions
#

"Returns index of range in which, or after, `k` is placed."
@inline searchsortedrangelast(rs::VecUnitRangesSortedSet{K}, k) where {K} = searchsortedlast(rs.rstarts, K(k); lt=<)
@inline searchsortedrangelast(rs::UnitRangesSortedSet{K}, k) where {K} = searchsortedlast(rs.ranges, K(k))
@inline searchsortedrangelast(rs::UnitRangesSortedSubSet0, k) = beforefirstindex(rs)
@inline searchsortedrangelast(rs::UnitRangesSortedSubSet1{K}, k) where {K} =
    K(k) >= rs.kstart ? rs.firstindex : beforefirstindex(rs)
@inline searchsortedrangelast(rs::UnitRangesSortedSubSet{K,TU,P,Tix}, k) where {K,TU,P<:VecUnitRangesSortedSet,Tix} =
    searchsortedlast(rs.parent.rstarts, K(k), rs.firstindex, rs.lastindex, Base.ord(isless, identity, false))
@inline function searchsortedrangelast(rs::UnitRangesSortedSubSet{K,TU,P,Tix}, k) where {K,TU,P<:UnitRangesSortedSet,Tix}
    ir_k = searchsortedlast(rs.parent.ranges, K(k))
    if compare(rs.parent.ranges, ir_k, beforefirstindex(rs)) == 1 && compare(rs.parent.ranges, ir_k, pastlastindex(rs)) == -1
        return ir_k
    else
        return beforefirstindex(rs)
    end
end

"Returns index of range in which, or before, `k` is placed."
@inline searchsortedrangefirst(rs::VecUnitRangesSortedSet{K}, k) where {K} = searchsortedfirst(rs.rstops, K(k); lt=isless)
@inline function searchsortedrangefirst(rs::UnitRangesSortedSet{K}, k) where {K}
    ir_k = searchsortedrangelast(rs, K(k))
    if ir_k != beforefirstindex(rs) && in(K(k), getindex(rs, ir_k))
        return ir_k
    else
        return advance(rs, ir_k)
    end
end
@inline searchsortedrangefirst(rs::UnitRangesSortedSubSet0, k) = pastlastindex(rs)
@inline searchsortedrangefirst(rs::UnitRangesSortedSubSet1{K}, k) where {K} = K(k) <= rs.kstop ? rs.lastindex : pastlastindex(rs)
@inline searchsortedrangefirst(rs::UnitRangesSortedSubSet{K,TU,P,Tix}, k) where {K,TU,P<:VecUnitRangesSortedSet,Tix} =
    searchsortedfirst(rs.parent.rstops, K(k), rs.firstindex, rs.lastindex, Base.ord(isless, identity, false))
@inline function searchsortedrangefirst(rs::UnitRangesSortedSubSet{K,TU,P,Tix}, k) where {K,TU,P<:UnitRangesSortedSet,Tix}
    ir_k = searchsortedrangelast(rs, K(k))
    if ir_k != beforefirstindex(rs) && in(K(k), getindex(rs, ir_k))
        return ir_k
    else
        return advance(rs, ir_k)
    end
end

"Returns range of `rs` indexes which coincide or concluded in `kk` range."
@inline searchsortedrange(rs::AbstractUnitRangesSortedContainer, kk::AbstractRange) =
    create_indexrange(rs, searchsortedrangefirst(rs, first(kk)), searchsortedrangelast(rs, last(kk)))
@inline searchsortedrange(rs::AbstractUnitRangesSortedSubSet, kk::AbstractRange) =
    create_indexrange(rs.parent, searchsortedrangefirst(rs, first(kk)), searchsortedrangelast(rs, last(kk)))

"Returns indexes of range in `rs` in which `k` may be inserted. Or negative range in the case of `k` is
 between `rs` ranges, and indices of resulted range is the indexes of that neighbors."
@inline function searchsortedrange(rs::AbstractUnitRangesSortedContainer{K}, k) where {K}
    ir_k = searchsortedrangelast(rs, K(k))
    if ir_k != beforefirstindex(rs) && in(K(k), getindex(rs, ir_k))
        return create_indexrange(rs, ir_k, ir_k)
    else
        return create_indexrange(rs, advance(rs, ir_k), ir_k)
    end
end


@inline index_status(rs::VecUnitRangesSortedSet, ir) = 0 <= ir <= length(rs)+1 ? 1 : 0
@inline index_status(rs::UnitRangesSortedSet, ir) = status((rs.ranges, ir))

@inline function Base.findfirst(pred::Function, rs::AbstractUnitRangesSortedSet)
    for i in eachindex(rs)
        pred(getindex(rs, i)) && return i
    end
    return nothing
end

Base.findall(pred::Function, rs::AbstractUnitRangesSortedSet) = collect(r for r in rs if pred(r))



#
#  Iterators
#

@inline Base.iterate(rs::UnitRangesSortedSubSet0) = nothing
@inline Base.iterate(rs::UnitRangesSortedSubSet0, state) = nothing
@inline Base.iterate(rrs::Base.Iterators.Reverse{T}) where {T<:UnitRangesSortedSubSet0} = nothing
@inline Base.iterate(rrs::Base.Iterators.Reverse{T}, state) where {T<:UnitRangesSortedSubSet0} = nothing


@inline Base.iterate(rs::UnitRangesSortedSubSet1) = (rs.singlerange, 1)
@inline Base.iterate(rs::UnitRangesSortedSubSet1, state) = nothing
@inline Base.iterate(rrs::Base.Iterators.Reverse{T}) where {T<:UnitRangesSortedSubSet1} =
    (rrs.itr.singlerange, 1)
@inline Base.iterate(rrs::Base.Iterators.Reverse{T}, state) where {T<:UnitRangesSortedSubSet1} = nothing


@inline Base.iterate(rs::AbstractUnitRangesSortedSubSet{K,TU,P}) where {K,TU,P<:VecUnitRangesSortedSet} =
    (rs.firstrange, rs.firstindex + 1)
@inline function Base.iterate(rs::AbstractUnitRangesSortedSubSet{K,TU,P}, state) where {K,TU,P<:VecUnitRangesSortedSet}
    if state < rs.lastindex
        return (getindex(rs.parent, state), state + 1)
    elseif state == rs.lastindex
        return (rs.lastrange, state + 1)
    else
        return nothing
    end
end
@inline function Base.iterate(rrs::Base.Iterators.Reverse{T}) where {T<:AbstractUnitRangesSortedSubSet{K,TU,P}} where
                                                                    {K,TU,P<:VecUnitRangesSortedSet}
    return (rrs.itr.lastrange, rrs.itr.lastindex - 1)
end
@inline function Base.iterate(rrs::Base.Iterators.Reverse{T}, state) where
                              {T<:AbstractUnitRangesSortedSubSet{K,TU,P}} where {K,TU,P<:VecUnitRangesSortedSet}
    if rrs.itr.firstindex < state
        return (getindex(rrs.itr.parent, state), state - 1)
    elseif rrs.itr.firstindex == state
        return (rrs.itr.firstrange, state - 1)
    else
        return nothing
    end
end


@inline Base.iterate(rs::AbstractUnitRangesSortedSubSet{K,TU,P}) where {K,TU,P<:UnitRangesSortedSet} =
    (rs.firstrange, advance((rs.parent.ranges, rs.firstindex)))
@inline function Base.iterate(rs::AbstractUnitRangesSortedSubSet{K,TU,P}, state) where {K,TU,P<:UnitRangesSortedSet}
    if compare(rs.parent.ranges, state, rs.lastindex) == -1
        return (getindex(rs.parent, state), advance((rs.parent.ranges, state)))
    elseif compare(rs.parent.ranges, state, rs.lastindex) == 0
        return (rs.lastrange, advance((rs.parent.ranges, state)))
    else
        return nothing
    end
end
@inline function Base.iterate(rrs::Base.Iterators.Reverse{T}) where {T<:AbstractUnitRangesSortedSubSet{K,TU,P}} where
                                                                    {K,TU,P<:UnitRangesSortedSet}
    return (rrs.itr.lastrange, regress((rrs.itr.parent.ranges, rrs.itr.lastindex)))
end
@inline function Base.iterate(rrs::Base.Iterators.Reverse{T}, state) where
                              {T<:AbstractUnitRangesSortedSubSet{K,TU,P}} where {K,TU,P<:UnitRangesSortedSet}
    if compare(rrs.itr.parent.ranges, rrs.itr.firstindex, state) == -1
        return (getindex(rrs.itr.parent, state), regress((rrs.itr.parent.ranges, state)))
    elseif compare(rrs.itr.parent.ranges, rrs.itr.firstindex, state) == 0
        return (rrs.itr.firstrange, regress((rrs.itr.parent.ranges, state)))
    else
        return nothing
    end
end



@inline function Base.iterate(rs::AbstractUnitRangesSortedSet, state = firstindex(rs))
    if state != pastlastindex(rs)
        return (getindex(rs, state), advance(rs, state))
    else
        return nothing
    end
end
@inline function Base.iterate(rrs::Base.Iterators.Reverse{T}, state = lastindex(rrs.itr)) where {T<:AbstractUnitRangesSortedSet}
    if state != beforefirstindex(rrs.itr)
        return (getindex(rrs.itr, state), regress(rrs.itr, state))
    else
        return nothing
    end
end


@inline function Base.iterate(ur::URSSIndexURange, state = (first(ur), 0))
    st, i = state
    if i < length(ur)
        return (st, (advance(ur.parent, st), i + 1))
    else
        return nothing
    end
end


@inline function eachindexiterate(rs::AbstractUnitRangesSortedSet, state = firstindex(rs))
    if state != pastlastindex(rs)
        return (state, advance(rs, state))
    else
        return nothing
    end
end
@inline function eachindexiterate(rrs::Base.Iterators.Reverse{T}, state = lastindex(rrs.itr)) where {T<:AbstractUnitRangesSortedSet}
    if state != beforefirstindex(rrs.itr)
        return (state, regress(rrs.itr, state))
    else
        return nothing
    end
end


struct EachIndexIterator{It}
    itr::It
end
"`eachindex` iterator"
@inline eachindexiterator(itr) = EachIndexIterator(itr)
@inline function Base.iterate(it::EachIndexIterator, state...)
    y = eachindexiterate(it.itr, state...)
    if y !== nothing
        return (y[1], y[2])
    else
        return nothing
    end
end
Base.eltype(::Type{EachIndexIterator{It}}) where {It} = eltype(It)
Base.IteratorEltype(::Type{EachIndexIterator{It}}) where {It} = Base.HasEltype()
Base.IteratorSize(::Type{<:EachIndexIterator}) = Base.HasShape{1}()
Base.length(it::EachIndexIterator) = length(it.itr)
Base.size(it::EachIndexIterator) = (length(it.itr),)
Iterators.reverse(it::EachIndexIterator) = EachIndexIterator(Iterators.reverse(it.itr))


Base.eachindex(rs::VecUnitRangesSortedSet) = eachindex(rs.rstarts)
Base.eachindex(rs::UnitRangesSortedSet) = eachindexiterator(rs)
Base.eachindex(rs::T) where {T<:AbstractUnitRangesSortedSubSet} = eachindexiterator(rs)
Base.eachindex(rrs::Base.Iterators.Reverse{T}) where {T<:AbstractUnitRangesSortedSet} = Iterators.reverse(eachindex(rrs.itr))
# TODO: check for https://github.com/JuliaLang/julia/pull/43110/files

#
# Assignments
#


@inline Base.in(kk::AbstractRange, rs::AbstractUnitRangesSortedSet{K,TU}) where {K,TU} = _in(to_urange(TU, kk), rs)
@inline Base.in(kk::TU, rs::AbstractUnitRangesSortedSet{K,TU}) where {K,TU<:AbstractRange} = _in(kk, rs)

@inline function _in(kk::TU, rs::AbstractUnitRangesSortedSet{K,TU}) where {K,TU}
    # fast check for cached range index
    if (ir = rs.lastusedrangeindex) != beforefirstindex(rs)
        if issubset(kk, getindex(rs, ir))
            return true
        end
    end
    # cached range' index miss (or index is not stored), thus try to search
    iir_kk = searchsortedrange(rs, kk)
    rs.lastusedrangeindex = last(iir_kk)
    if length(iir_kk) != 1
        return false
    elseif issubset(kk, getindex(rs, first(iir_kk)))
        return true
    else
        return false
    end
end

@inline function Base.in(key, rs::AbstractUnitRangesSortedSet{K}) where K
    k = K(key)
    # fast check for cached range index
    if (ir = rs.lastusedrangeindex) != beforefirstindex(rs)
        r_start, r_stop = getindex_tuple(rs, ir)
        if r_start <= k <= r_stop
            return true
        end
    end
    # cached range index miss (or index not stored), thus try search
    ir_k = searchsortedrangelast(rs, k)
    if ir_k != beforefirstindex(rs)  # `k` is not before the start of first range
        r_start, r_stop = getindex_tuple(rs, ir_k)
        if k <= r_stop  # is `k` inside of range
            rs.lastusedrangeindex = ir_k
            return true
        end
    end
    rs.lastusedrangeindex = beforefirstindex(rs)
    return false
end

@inline Base.in(kk::AbstractRange, su::UnitRangesSortedSubSet0) = false
@inline Base.in(key, su::UnitRangesSortedSubSet0) = false
@inline Base.in(kk::AbstractRange, su::UnitRangesSortedSubSet1{K,TU}) where {K,TU} = issubset(to_urange(TU, kk), su.singlerange)
@inline Base.in(kk::TU, su::UnitRangesSortedSubSet1{K,TU}) where {K,TU<:AbstractRange} = issubset(kk, su.singlerange)
@inline Base.in(key, su::UnitRangesSortedSubSet1{K}) where {K} = in(K(key), su.singlerange)

@inline Base.in(kk::AbstractRange, su::UnitRangesSortedSubSet{K,TU}) where {K,TU} = _in(to_urange(TU, kk), su)
@inline Base.in(kk::TU, su::UnitRangesSortedSubSet{K,TU}) where {K,TU<:AbstractRange} = _in(kk, su)
@inline function _in(kk::TU, su::UnitRangesSortedSubSet{K,TU,P}) where {K,TU,P<:VecUnitRangesSortedSet}
    # bounds check
    (su.kstart <= first(kk) <= last(kk) <= su.kstop) || return false

    # fast check for cached range index
    if (ir = su.parent.lastusedrangeindex) != beforefirstindex(su.parent) &&
        su.firstindex <= first(ir) <= last(ir) <= su.lastindex
        if issubset(kk, getindex(su, ir))
            return true
        end
    end
    # cached range' index miss (or index is not stored), thus try search
    iir_kk = searchsortedrange(su, kk)
    su.parent.lastusedrangeindex = last(iir_kk)
    if length(iir_kk) != 1
        return false
    elseif issubset(kk, getindex(su, first(iir_kk)))
        return true
    else
        return false
    end
end
@inline function _in(kk::TU, su::UnitRangesSortedSubSet{K,TU}) where {K,TU}
    # bounds check
    issubset(kk, to_urange(TU, su.kstart, su.kstop)) || return false

    iir_kk = searchsortedrange(su, kk)
    su.parent.lastusedrangeindex = last(iir_kk)
    if length(iir_kk) != 1
        return false
    elseif issubset(kk, getindex(su, first(iir_kk)))
        return true
    else
        return false
    end
end

@inline function Base.in(key, su::UnitRangesSortedSubSet{K}) where {K}
    k = K(key)
    # check outbound
    (k < su.kstart || su.kstop < k) && return false

    # fast check for cached range index
    if (ir = su.parent.lastusedrangeindex) != beforefirstindex(su.parent) &&
        indexcompare(su, su.firstindex, ir) < 1 &&
        indexcompare(su, ir, su.lastindex) < 1
        r_start, r_stop = getindex_tuple(su, ir)
        if r_start <= k <= r_stop
            return true
        end
    end
    # cached range index miss (or index not stored), thus try search
    ir_k = searchsortedrangelast(su, k)
    if ir_k != beforefirstindex(su)  # `k` is not before the start of first range
        r_start, r_stop = getindex_tuple(su, ir_k)
        if r_start <= k <= r_stop  # is `k` inside of range
            su.parent.lastusedrangeindex = ir_k
            return true
        end
    end
    su.parent.lastusedrangeindex = beforefirstindex(su.parent)
    return false
end



function Base.push!(rs::VecUnitRangesSortedSet{K,TU}, II::Union{AbstractVector,AbstractSet,NTuple}) where {K,TU}
    for r in II
        push!(rs, r)
    end
    rs
end


@inline Base.push!(rs::VecUnitRangesSortedSet{K,TU}, kk::AbstractRange) where {K,TU} = _push!(rs, to_urange(TU, kk))
@inline Base.push!(rs::VecUnitRangesSortedSet{K,TU}, kk::TU) where {K,TU<:AbstractRange} = _push!(rs, kk)

function _push!(rs::VecUnitRangesSortedSet{K,TU}, kk::TU) where {K,TU}

    if length(kk) == 0
        return rs
    elseif length(kk) == 1
        push!(rs, first(kk))
        return rs
    end

    rs.lastusedrangeindex = beforefirstindex(rs)

    iir_kk = searchsortedrange(rs, kk)

    # `kk` already exist in `rs` and the same as the corresponding range in `rs`, nothing to do
    if length(iir_kk) == 1 && kk == getindex(rs, first(iir_kk))
        return rs

    # delete all inside `kk` ranged in `rs`
    else
        delete!(rs, kk)
    end

    # get neighbors pointers
    # Note: the range will be zero-length, thus
    # `first(iir_kk)` will be point to range in `rs` on right side to `kk` and
    # `last(iir_kk)` will be point to range in `rs` on left side to `kk`.
    iir_kk = searchsortedrange(rs, first(kk))
    @boundscheck @assert iir_kk == searchsortedrange(rs, last(kk)) "FIXME: Something went wrong."

    iir_kk_left = last(iir_kk)
    iir_kk_right = first(iir_kk)

    # `kk` is adjoined in both sides with `rs` ranges, thus join them all
    if iir_kk_left != beforefirstindex(rs) && getindex_rangestop(rs, iir_kk_left) + 1 == first(kk) &&
       iir_kk_right != pastlastindex(rs) && getindex_rangestart(rs, iir_kk_right) - 1 == last(kk)
        rs.rstops[iir_kk_left] = getindex_rangestop(rs, iir_kk_right)
        deleteat!(rs.rstarts, iir_kk_right)
        deleteat!(rs.rstops, iir_kk_right)

    # `kk` is adjoin with `rs` range on left side, thus append to left range
    elseif iir_kk_left != beforefirstindex(rs) && getindex_rangestop(rs, iir_kk_left) + 1 == first(kk) &&
           iir_kk_right != pastlastindex(rs) && getindex_rangestart(rs, iir_kk_right) - 1 != last(kk)
        rs.rstops[iir_kk_left] = last(kk)

    # `kk` is adjoin with `rs` range on right side, thus prepend to right range
    elseif iir_kk_left != beforefirstindex(rs) && getindex_rangestop(rs, iir_kk_left) + 1 != first(kk) &&
           iir_kk_right != pastlastindex(rs) && getindex_rangestart(rs, iir_kk_right) - 1 == last(kk)
        rs.rstarts[iir_kk_right] = first(kk)

    # `kk` is separate from both sides, insert it
    elseif iir_kk_left != beforefirstindex(rs) && getindex_rangestop(rs, iir_kk_left) + 1 != first(kk) &&
           iir_kk_right != pastlastindex(rs) && getindex_rangestart(rs, iir_kk_right) - 1 != last(kk)
        insert!(rs.rstarts, iir_kk_right, first(kk))
        insert!(rs.rstops, iir_kk_right, last(kk))

    # `kk` is first range in `rs` and adjoin with range on right side, thus prepend `kk` to right range
    elseif iir_kk_left == beforefirstindex(rs) &&
           iir_kk_right != pastlastindex(rs) && getindex_rangestart(rs, iir_kk_right) - 1 == last(kk)
        rs.rstarts[iir_kk_right] = first(kk)

    # `kk` is separate first range in `rs`, insert it
    elseif iir_kk_left == beforefirstindex(rs) &&
           iir_kk_right != pastlastindex(rs) && getindex_rangestart(rs, iir_kk_right) - 1 != last(kk)
        insert!(rs.rstarts, iir_kk_right, first(kk))
        insert!(rs.rstops, iir_kk_right, last(kk))

    # `kk` is last range in `rs` and adjoin with range on left side, thus append `kk` to left range
    elseif iir_kk_left != beforefirstindex(rs) && getindex_rangestop(rs, iir_kk_left) + 1 == first(kk) &&
           iir_kk_right == pastlastindex(rs)
        rs.rstops[iir_kk_left] = last(kk)

    # `kk` is separate last range in `rs`, insert it
    elseif iir_kk_left != beforefirstindex(rs) && getindex_rangestop(rs, iir_kk_left) + 1 != first(kk) &&
           iir_kk_right == pastlastindex(rs)
        insert!(rs.rstarts, iir_kk_right, first(kk))
        insert!(rs.rstops, iir_kk_right, last(kk))

    # `rs` is empty
    elseif iir_kk_left == beforefirstindex(rs) && iir_kk_right == pastlastindex(rs)
        insert!(rs.rstarts, 1, first(kk))
        insert!(rs.rstops, 1, last(kk))

    else
        throw(AssertionError("FIXME: Should be unreachable."))
    end

    return rs
end

function Base.push!(rs::UnitRangesSortedSet{K,TU}, II::Union{AbstractVector,AbstractSet,NTuple}) where {K,TU}
    for r in II
        push!(rs, r)
    end
    rs
end

@inline Base.push!(rs::UnitRangesSortedSet{K,TU}, kk::AbstractRange) where {K,TU} = _push!(rs, to_urange(TU, kk))
@inline Base.push!(rs::UnitRangesSortedSet{K,TU}, kk::TU) where {K,TU<:AbstractRange} = _push!(rs, kk)
function _push!(rs::UnitRangesSortedSet{K,TU}, kk::TU) where {K,TU}

    if length(kk) == 0
        return rs
    elseif length(kk) == 1
        push!(rs, first(kk))
        return rs
    end

    rs.lastusedrangeindex = beforefirstindex(rs)

    iir_kk = searchsortedrange(rs, kk)

    # `kk` already exist in `rs` and the same as the corresponding range in `rs`, nothing to do
    if length(iir_kk) == 1 && kk == getindex(rs, first(iir_kk))
        return rs

    # delete all inside `kk` ranged in `rs`
    else
        delete!(rs, kk)
    end

    # get neighbors pointers
    # Note: the range will be zero-length, thus
    # `first(iir_kk)` will be point to range in `rs` on right side to `kk` and
    # `last(iir_kk)` will be point to range in `rs` on left side to `kk`.
    iir_kk = searchsortedrange(rs, first(kk))
    @boundscheck @assert iir_kk == searchsortedrange(rs, last(kk)) "FIXME: Something went wrong."

    iir_kk_left = last(iir_kk)
    iir_kk_right = first(iir_kk)
    iir_kk_left_rangestop = iir_kk_left != beforefirstindex(rs) ? getindex_rangestop(rs, iir_kk_left) : K(0)
    iir_kk_right_rangestart = iir_kk_right != pastlastindex(rs) ? getindex_rangestart(rs, iir_kk_right) : K(0)
    #iir_kk_left_rangestop = iir_kk_left != beforefirstindex(rs) ?
    #                                             getindex_rangestop(rs, iir_kk_left) : typemin(K)
    #iir_kk_right_rangestart = iir_kk_right != pastlastindex(rs) ?
    #                                                getindex_rangestart(rs, iir_kk_right) : typemax(K)

    # `kk` is adjoined in both sides with `rs` ranges, thus join them all
    if iir_kk_left != beforefirstindex(rs) && iir_kk_left_rangestop + 1 == first(kk) &&
       iir_kk_right != pastlastindex(rs) && iir_kk_right_rangestart - 1 == last(kk)
        rs.ranges[iir_kk_left] = getindex_rangestop(rs, iir_kk_right)
        delete!((rs.ranges, iir_kk_right))

    # `kk` is adjoin with `rs` range on left side, thus append to left range
    elseif iir_kk_left != beforefirstindex(rs) && iir_kk_left_rangestop + 1 == first(kk) &&
           iir_kk_right != pastlastindex(rs) && iir_kk_right_rangestart - 1 != last(kk)
        rs.ranges[iir_kk_left] = last(kk)

    # `kk` is adjoin with `rs` range on right side, thus prepend to right range
    elseif iir_kk_left != beforefirstindex(rs) && iir_kk_left_rangestop + 1 != first(kk) &&
           iir_kk_right != pastlastindex(rs) && iir_kk_right_rangestart - 1 == last(kk)
        rs.ranges[first(kk)] = getindex_rangestop(rs, iir_kk_right)
        delete!((rs.ranges, iir_kk_right))

    # `kk` is separate from both sides, insert it
    elseif iir_kk_left != beforefirstindex(rs) && iir_kk_left_rangestop + 1 != first(kk) &&
           iir_kk_right != pastlastindex(rs) && iir_kk_right_rangestart - 1 != last(kk)
        rs.ranges[first(kk)] = last(kk)

    # `kk` is first range in `rs` and adjoin with range on right side, thus prepend `kk` to right range
    elseif iir_kk_left == beforefirstindex(rs) &&
           iir_kk_right != pastlastindex(rs) && iir_kk_right_rangestart - 1 == last(kk)
        rs.ranges[first(kk)] = getindex_rangestop(rs, iir_kk_right)
        delete!((rs.ranges, iir_kk_right))

    # `kk` is separate first range in `rs`, insert it
    elseif iir_kk_left == beforefirstindex(rs) &&
           iir_kk_right != pastlastindex(rs) && iir_kk_right_rangestart - 1 != last(kk)
        rs.ranges[first(kk)] = last(kk)

    # `kk` is last range in `rs` and adjoin with range on left side, thus append `kk` to left range
    elseif iir_kk_left != beforefirstindex(rs) && iir_kk_left_rangestop + 1 == first(kk) &&
           iir_kk_right == pastlastindex(rs)
        rs.ranges[iir_kk_left] = last(kk)

    # `kk` is separate last range in `rs`, insert it
    elseif iir_kk_left != beforefirstindex(rs) && iir_kk_left_rangestop + 1 != first(kk) &&
           iir_kk_right == pastlastindex(rs)
        rs.ranges[first(kk)] = last(kk)

    # `rs` is empty
    elseif iir_kk_left == beforefirstindex(rs) && iir_kk_right == pastlastindex(rs)
        rs.ranges[first(kk)] = last(kk)

    else
        throw(AssertionError("FIXME: Should be unreachable."))
    end

    return rs
end

@inline function check_exist_and_update(rs, k)

    # fast check for cached range index
    if (ir = rs.lastusedrangeindex) != beforefirstindex(rs)
        r_start, r_stop = getindex_tuple(rs, ir)
        if r_start <= k <= r_stop
            return nothing
        end
    end

    ir_k = searchsortedrangelast(rs, k)
    #
    #sstatus = status((rs.ranges, ir_k))
    @boundscheck if index_status(rs, ir_k) == 0 # invalid range index
        throw(BoundsError(rs, k))
    end

    # check the index exist and update its data
    if ir_k != beforefirstindex(rs)  # `k` is not before the first index
        if k <= getindex_rangestop(rs, ir_k)
            rs.lastusedrangeindex = ir_k
            return nothing
        end
    end

    return ir_k
end


function Base.push!(rs::VecUnitRangesSortedSet{K,TU}, key) where {K,TU}
    k = K(key)

    ir_k = check_exist_and_update(rs, k)
    ir_k == nothing && return rs

    if length(rs) == 0
        push!(rs.rstarts, k)
        push!(rs.rstops, k)
        rs.lastusedrangeindex = 1
        return rs
    end

    if ir_k == beforefirstindex(rs)  # `k` is before the first range
        if rs.rstarts[1] - k > 1  # there is will be gap in indices after inserting
            pushfirst!(rs.rstarts, k)
            pushfirst!(rs.rstops, k)
        else  # prepend to first range
            rs.rstarts[1] = k
        end
        rs.lastusedrangeindex = 1
        return rs
    end

    r = getindex(rs, ir_k)

    if k >= rs.rstarts[end]  # `k` is after the last range start
        if k > last(r) + 1  # there is will be the gap in indices after inserting
            push!(rs.rstarts, k)
            push!(rs.rstops, k)
        else  # just append to last range
            rs.rstops[ir_k] += 1
        end
        rs.lastusedrangeindex = lastindex(rs)
        return rs
    end

    # `k` is somewhere between indices
    stnext = advance(rs, ir_k)
    rnext = getindex(rs, stnext)

    if first(rnext) - last(r) == 2  # join ranges
        rs.rstarts[ir_k] = first(r)
        rs.rstops[ir_k] = last(rnext)
        deleteat!(rs.rstarts, stnext)
        deleteat!(rs.rstops, stnext)
        rs.lastusedrangeindex = ir_k
    elseif k - last(r) == 1  # append to left range
        rs.rstops[ir_k] += 1
        rs.lastusedrangeindex = ir_k
    elseif first(rnext) - k == 1  # prepend to right range
        rs.rstarts[stnext] -= 1
        rs.lastusedrangeindex = stnext
    else  # insert single element range
        insert!(rs.rstarts, stnext, k)
        insert!(rs.rstops, stnext, k)
        rs.lastusedrangeindex = stnext
    end

    return rs

end


function Base.push!(rs::UnitRangesSortedSet{K,TU}, key) where {K,TU}
    k = K(key)

    ir_k = check_exist_and_update(rs, k)
    ir_k == nothing && return rs

    if length(rs) == 0
        rs.ranges[k] = k
        rs.lastusedrangeindex = firstindex(rs)
        return rs
    end

    if ir_k == beforefirstindex(rs)  # `k` is before the first index
        stnext = advance(rs, ir_k)
        rnext = getindex(rs, stnext)
        if first(rnext) - k > 1  # there is will be gap in indices after inserting
            rs.ranges[k] = k
        else  # prepend to first range
            rs.ranges[k] = last(rnext)
            delete!((rs.ranges, stnext))
        end
        rs.lastusedrangeindex = firstindex(rs)
        return rs
    end

    r = getindex(rs, ir_k)

    if k > last(last(rs)) # `k` is after the last range end index
        if last(r) + 1 < k  # there is will be the gap in indices after inserting
            rs.ranges[k] = k
        else  # just append to last range
            rs.ranges[ir_k] = last(r)+1
        end
        rs.lastusedrangeindex = lastindex(rs)
        return rs
    end

    rs.lastusedrangeindex = beforefirstindex(rs)

    # `k` is somewhere between indices
    stnext = advance(rs, ir_k)
    rnext = getindex(rs, stnext)

    if first(rnext) - last(r) == 2  # join ranges
        rs.ranges[ir_k] = last(rnext)
        delete!((rs.ranges, stnext))
    elseif k - last(r) == 1  # append to left range
        rs.ranges[ir_k] = k
        rs.lastusedrangeindex = ir_k
    elseif first(rnext) - k == 1  # prepend to right range
        rs.ranges[k] = last(rnext)
        delete!((rs.ranges, stnext))
    else  # insert single element range
        rs.ranges[k] = k
    end

    return rs

end


function Base.delete!(rs::VecUnitRangesSortedSet{K,TU}, II::T) where {K,TU, T<:Union{AbstractVector,AbstractSet,NTuple}}
    if hasmethod(Iterators.reverse, Tuple{T}) &&
       hasmethod(keys, Tuple{T})
        for r in Iterators.reverse(II)
            delete!(rs, r)
        end
    else
        for r in II
            delete!(rs, r)
        end
    end
    rs
end

@inline Base.delete!(rs::VecUnitRangesSortedSet{K,TU}, kk::AbstractRange) where {K,TU} = _delete!(rs, to_urange(TU, kk))
@inline Base.delete!(rs::VecUnitRangesSortedSet{K,TU}, kk::TU) where {K,TU<:AbstractRange} = _delete!(rs, kk)
function _delete!(rs::VecUnitRangesSortedSet{K,TU}, kk::TU) where {K,TU}
    length(kk) == 0 && return rs

    iir_kk = searchsortedrange(rs, kk)
    length(iir_kk) == 0 && return rs

    rs.lastusedrangeindex = beforefirstindex(rs)

    # delete all concluded ranges
    todelete1 = getindex_rangestart(rs, first(iir_kk)) < first(kk) ? first(iir_kk) + 1 : first(iir_kk)
    todelete2 = getindex_rangestop(rs, last(iir_kk)) > last(kk) ? last(iir_kk) - 1 : last(iir_kk)
    splice!(rs.rstarts, todelete1:todelete2)
    splice!(rs.rstops, todelete1:todelete2)

    # remains the side ranges to delete
    # Note: there is not possible `beforefirstindex` or `pastlastindex` in `iir_kk`
    iir_kk = searchsortedrange(rs, kk)
    length(iir_kk) == 0 && return rs

    iir_kk_rangestart = getindex_rangestart(rs, first(iir_kk))
    iir_kk_rangestop = getindex_rangestop(rs, last(iir_kk))

    # `kk` intersects only one range in `rs`
    if length(iir_kk) == 1

        # inside one range, thus split it
        if iir_kk_rangestart < first(kk) && last(kk) < iir_kk_rangestop
            insert!(rs.rstarts, advance(rs, first(iir_kk)), last(kk) + 1)
            insert!(rs.rstops, first(iir_kk), first(kk) - 1)

        # `kk` intersects with, or inside in, one range from left side, thus shrink from left side
        elseif iir_kk_rangestart >= first(kk) && last(kk) < iir_kk_rangestop
            rs.rstarts[first(iir_kk)] = last(kk) + 1

        # `kk` intersects with, or inside in, one range from right side, thus shrink from right side
        # inside one range, ended to left side, thus shrink from right side
        elseif iir_kk_rangestart < first(kk) && last(kk) >= iir_kk_rangestop
            rs.rstops[first(iir_kk)] = first(kk) - 1

        else
            throw(AssertionError("FIXME: Should be unreachable."))
        end

    # All remaining cases are with two ranges from `rs`,
    # and both side ranges only partly intersects with `kk`.
    elseif length(iir_kk) == 2

        # shrink second (right) range from the left side
        rs.rstarts[last(iir_kk)] = last(kk) + 1

        # shrink first (left) range from the right side
        rs.rstops[first(iir_kk)] = first(kk) - 1

    else
        throw(AssertionError("FIXME: Should be unreachable."))
    end

    return rs
end


function Base.delete!(rs::UnitRangesSortedSet{K,TU}, II::Union{AbstractVector,AbstractSet,NTuple})  where {K,TU}
    for r in II
        delete!(rs, r)
    end
    rs
end

@inline Base.delete!(rs::UnitRangesSortedSet{K,TU}, kk::AbstractRange) where {K,TU} = _delete!(rs, to_urange(TU, kk))
@inline Base.delete!(rs::UnitRangesSortedSet{K,TU}, kk::TU) where {K,TU<:AbstractRange} = _delete!(rs, kk)
function _delete!(rs::UnitRangesSortedSet{K,TU}, kk::TU) where {K,TU}
    length(kk) == 0 && return rs

    iir_kk = searchsortedrange(rs, kk)
    length(iir_kk) == 0 && return rs

    rs.lastusedrangeindex = beforefirstindex(rs)

    # delete all concluded ranges
    todelete1 = getindex_rangestart(rs, first(iir_kk)) < first(kk) ? advance(rs, first(iir_kk)) : first(iir_kk)
    todelete2 = getindex_rangestop(rs, last(iir_kk)) > last(kk) ? regress(rs, last(iir_kk)) : last(iir_kk)
    for st in onlysemitokens(inclusive(rs.ranges, todelete1, todelete2))
        delete!((rs.ranges, st))
    end

    # remains the side ranges to delete
    # Note: there is not possible `beforefirstindex` or `pastlastindex` in `iir_kk`
    iir_kk = searchsortedrange(rs, kk)
    length(iir_kk) == 0 && return rs

    iir_kk_rangestart = getindex_rangestart(rs, first(iir_kk))
    iir_kk_rangestop = getindex_rangestop(rs, last(iir_kk))

    # `kk` intersects only one range in `rs`
    if length(iir_kk) == 1

        # inside one range, thus split it
        if iir_kk_rangestart < first(kk) && last(kk) < iir_kk_rangestop
            rs.ranges[last(kk) + 1] = iir_kk_rangestop
            rs.ranges[first(iir_kk)] = first(kk) - 1

        # `kk` intersects with, or inside in, one range from left side, thus shrink from left side
        elseif iir_kk_rangestart >= first(kk) && last(kk) < iir_kk_rangestop
            rs.ranges[last(kk) + 1] = iir_kk_rangestop
            delete!((rs.ranges, first(iir_kk)))

        # `kk` intersects with, or inside in, one range from right side, thus shrink from right side
        # inside one range, ended to left side, thus shrink from right side
        elseif iir_kk_rangestart < first(kk) && last(kk) >= iir_kk_rangestop
            rs.ranges[first(iir_kk)] = first(kk) - 1

        else
            throw(AssertionError("FIXME: Should be unreachable."))
        end

    # All remaining cases are with two ranges from `rs`,
    # and both side ranges only partly intersects with `kk`.
    elseif length(iir_kk) == 2

        # shrink second (right) range from the left side
        rs.ranges[last(kk) + 1] = iir_kk_rangestop
        delete!((rs.ranges, last(iir_kk)))


        # shrink first (left) range from the right side
        rs.ranges[first(iir_kk)] = first(kk) - 1

    else
        throw(AssertionError("FIXME: Should be unreachable."))
    end

    return rs
end


function Base.delete!(rs::VecUnitRangesSortedSet{K,TU}, key) where {K,TU}
    k = K(key)

    length(rs) == 0 && return rs

    ir_k = searchsortedrangelast(rs, k)

    if ir_k == beforefirstindex(rs)  # `k` is before first range
        return rs
    end

    r = getindex(rs, ir_k)

    if k > last(r)  # `k` is outside of range
        return rs
    end

    if last(r) - first(r) + 1 == 1
        deleteat!(rs.rstarts, ir_k)
        deleteat!(rs.rstops, ir_k)
    elseif k == last(r)  # last index in range
        rs.rstops[ir_k] -= 1
    elseif k == first(r)  # first index in range
        rs.rstarts[ir_k] += 1
    else
        insert!(rs.rstarts, advance(rs, ir_k), k+1)
        insert!(rs.rstops, advance(rs, ir_k), last(r))
        rs.rstops[ir_k] = k-1
    end

    rs.lastusedrangeindex = beforefirstindex(rs)

    return rs
end


function Base.delete!(rs::UnitRangesSortedSet{K,TU}, key) where {K,TU}
    k = K(key)

    length(rs) == 0 && return rs

    ir_k = searchsortedrangelast(rs, k)

    if ir_k == beforefirstindex(rs)  # `k` is before first index
        return rs
    end

    r = getindex(rs, ir_k)

    if k > last(r)  # `k` is outside of range
        return rs
    end

    if last(r) - first(r) + 1 == 1
        delete!((rs.ranges, ir_k))
    elseif k == last(r)  # last index in range
        rs.ranges[ir_k] = last(r) - 1
    elseif k == first(r)  # first index in range
        delete!((rs.ranges, ir_k))
        rs.ranges[k+1] = last(r)
    else
        rs.ranges[ir_k] = k - 1
        rs.ranges[k+1] = last(r)
    end

    rs.lastusedrangeindex = beforefirstindex(rs)

    return rs
end


function Base.pop!(rs::AbstractUnitRangesSortedContainer, k)
    if (r = getrange(rs, k)) !== nothing
        delete!(rs, k)
        return r
    else
        throw(KeyError(k))
    end
end
function Base.pop!(rs::AbstractUnitRangesSortedContainer, k, default)
    if (r = getrange(rs, k)) !== nothing
        delete!(rs, k)
        return r
    else
        return default
    end
end


Base.empty(rs::T) where {T<:AbstractUnitRangesSortedContainer} = T()

function Base.empty!(rs::VecUnitRangesSortedSet)
    empty!(rs.rstarts)
    empty!(rs.rstops)
    rs.lastusedrangeindex = beforefirstindex(rs)
    rs
end
function Base.empty!(rs::AbstractUnitRangesSortedContainer)
    empty!(rs.ranges)
    rs.lastusedrangeindex = beforefirstindex(rs)
    rs
end

@inline Base.copy(rs::T) where {T<:VecUnitRangesSortedSet} = T(rs.lastusedrangeindex, copy(rs.rstarts), copy(rs.rstops))
@inline Base.copy(rs::T) where {T<:UnitRangesSortedSet} = T(rs.lastusedrangeindex, packcopy(rs.ranges))



@inline Base.union(rs::AbstractUnitRangesSortedSet, rss...) = union!(copy(rs), rss...)
@inline Base.union(rs::AbstractUnitRangesSortedSet, rs2) = union!(copy(rs), rs2)

@inline Base.union!(rs::AbstractUnitRangesSortedContainer, rss...) = union!(union!(rs, rss[1]), Base.tail(rss)...)
function Base.union!(rs::AbstractUnitRangesSortedContainer, rs2)
    issubset(rs2, rs) && return rs
    for r in rs2
        push!(rs, r)
    end
    return rs
end

@inline Base.union!(rs::AbstractUnitRangesSortedContainer, r::AbstractRange) = push!(rs, r)


@inline Base.setdiff(rs::AbstractUnitRangesSortedSet, rss...) = setdiff!(copy(rs), rss...)
@inline Base.setdiff(rs::AbstractUnitRangesSortedSet, rs2) = setdiff!(copy(rs), rs2)

@inline Base.setdiff!(rs::AbstractUnitRangesSortedContainer, rss...) = setdiff!(setdiff!(rs, rss[1]), Base.tail(rss)...)
@inline Base.setdiff!(rs::AbstractUnitRangesSortedContainer, rs2::AbstractRange) = delete!(rs, rs2)
function Base.setdiff!(rs::AbstractUnitRangesSortedContainer, rs2)
    if hasmethod(Iterators.reverse, Tuple{typeof(rs2)}) &&
       hasmethod(keys, Tuple{typeof(rs2)})
        for r in Iterators.reverse(rs2)
            delete!(rs, r)
        end
    else
        for r in rs2
            delete!(rs, r)
        end
    end
    return rs
end

@inline Base.symdiff(rs::AbstractUnitRangesSortedSet, sets...) = symdiff!(copy(rs), sets...)
@inline Base.symdiff(rs::AbstractUnitRangesSortedSet, s) = symdiff!(copy(rs), s)

@inline Base.symdiff!(rs::AbstractUnitRangesSortedContainer, rss...) = symdiff!(_symdiff!(rs, rss[1]), Base.tail(rss)...)
@inline Base.symdiff!(rs1::AbstractUnitRangesSortedContainer, rs2) = _symdiff!(rs1, rs2)
@inline Base.symdiff!(rs1::AbstractUnitRangesSortedContainer, rs2::AbstractSet) = _symdiff!(rs1, rs2) # AbstractSet ambiguity
function _symdiff!(rs1::AbstractUnitRangesSortedContainer{K,TU}, rs2) where {K,TU}
    for r in rs2
        ss = subset(rs1, r)
        vv = collect(ss)
        if length(vv) == 0
            push!(rs1, r)
        elseif length(vv) == 1
            delete!(rs1, vv[1])
            push!(rs1, to_urange(TU, first(r), safe_sub(first(vv[1]), 1)))
            push!(rs1, to_urange(TU, safe_add(last(vv[1]), 1), last(r)))
        else
            push!(rs1, to_urange(TU, first(r), safe_sub(first(vv[1]), 1)))
            push!(rs1, to_urange(TU, safe_add(last(vv[end]), 1), last(r)))
            rnext, rhead = Iterators.peel(Iterators.reverse(vv))
            delete!(rs1, rnext)
            for rr in rhead
                delete!(rs1, rr)
                push!(rs1, to_urange(TU, last(rr)+1, first(rnext)-1))
                rnext = rr
            end
        end
    end
    return rs1
end


@inline Base.:(==)(rs1::T1, rs2::T2) where {T1<:AbstractUnitRangesSortedSet, T2<:AbstractUnitRangesSortedSet} =
    T1 == T2 && isequal(rs1, rs2)

function Base.isequal(rs1::AbstractUnitRangesSortedSet, rs2::AbstractUnitRangesSortedSet)
    length(rs1) == length(rs2) || return false
    for (r1,r2) in zip(rs1, rs2)
        r1 == r2 || return false
    end
    return true
end

function Base.issubset(rs1::AbstractUnitRangesSortedSet, rs2::AbstractUnitRangesSortedSet)
    for r1 in rs1
        in(r1, rs2) || return false
    end
    return true
end
function Base.issubset(rs1::Union{AbstractSet,AbstractVector,AbstractRange,Tuple}, rs2::AbstractUnitRangesSortedSet)
    for r1 in rs1
        in(r1, rs2) || return false
    end
    return true
end
@inline Base.issubset(rs1::AbstractUnitRangesSortedSet, rs2::Union{AbstractSet,AbstractVector,AbstractRange,Tuple}) = _issubset1(rs1, rs2)
@inline Base.issubset(rs1::AbstractUnitRangesSortedSet, rs2::SortedSet) = _issubset1(rs1, rs2)  # ambiguity
function _issubset1(rs1::AbstractUnitRangesSortedSet, rs2)
    for r1 in rs1
        issubset(r1, rs2) || return false
    end
    return true
end
@inline Base.issubset(rs1::AbstractUnitRangesSortedSet, rs2::AbstractSet{T}) where {T<:AbstractRange} = _issubset2(rs1, rs2)
@inline Base.issubset(rs1::AbstractUnitRangesSortedSet, rs2::SortedSet{T}) where {T<:AbstractRange} = _issubset2(rs1, rs2)  # ambiguity
@inline Base.issubset(rs1::AbstractUnitRangesSortedSet, rs2::AbstractVector{T}) where {T<:AbstractRange} = _issubset2(rs1, rs2)
@inline Base.issubset(rs1::AbstractUnitRangesSortedSet, rs2::Tuple{Vararg{AbstractRange}}) = _issubset2(rs1, rs2)
function _issubset2(rs1::AbstractUnitRangesSortedSet, rs2)
    for r1 in rs1
        findfirst(s->issubset(r1, s), rs2) !== nothing || return false
    end
    return true
end


function Base.filter(pred::Function, rs::AbstractUnitRangesSortedContainer)
    res = empty(rs)
    for r in rs
        pred(r) && push!(res, r)
    end
    res
end
function Base.filter!(pred::Function, rs::AbstractUnitRangesSortedContainer)
    for r in Iterators.reverse(rs)
        pred(r) || delete!(rs, r)
    end
    rs
end


@inline Base.intersect!(rs::AbstractUnitRangesSortedContainer{K,TU}, kk::AbstractRange) where {K,TU} = __intersect!(rs, to_urange(TU,kk))
@inline Base.intersect!(rs::AbstractUnitRangesSortedContainer{K,TU}, kk::TU) where {K,TU<:AbstractRange} = __intersect!(rs, kk)
function __intersect!(rs::AbstractUnitRangesSortedContainer{K,TU}, kk::TU) where {K,TU}
    length(kk) == 0 && return empty!(rs)
    delete!(rs, safe_add(last(kk), 1):getindex_rangestop(rs, lastindex(rs)))
    delete!(rs, getindex_rangestart(rs, firstindex(rs)):safe_sub(first(kk), 1))
    rs
end

function Base.intersect!(rs1::AbstractUnitRangesSortedContainer, rs2::AbstractUnitRangesSortedSet)
    length(rs2) == 0 && return empty!(rs1)

    # cut both sides up to `rs2` indices
    intersect!(rs1, getindex_rangestart(rs2, firstindex(rs2)):getindex_rangestop(rs2, lastindex(rs2)))

    length(rs2) == 1 && return rs1

    if hasmethod(Iterators.reverse, Tuple{typeof(rs2)}) &&
       hasmethod(keys, Tuple{typeof(rs2)})
        rnext, rs2head = Iterators.peel(Iterators.reverse(rs2))
        for r in rs2head
            delete!(rs1, last(r)+1:first(rnext)-1)
            rnext = r
        end
    else
        rprev, rs2tail = Iterators.peel(rs2)
        for r in rs2tail
            delete!(rs1, last(rprev)+1:first(r)-1)
            rprev = r
        end
    end

    rs1
end

# AbstractSet ambiguity resolving
Base.intersect!(rs1::AbstractUnitRangesSortedContainer, rs2::AbstractSet) = _intersect1!(rs1, rs2)
Base.intersect!(rs1::AbstractUnitRangesSortedContainer, rs2::AbstractVector) = _intersect1!(rs1, rs2)
function _intersect1!(rs1::AbstractUnitRangesSortedContainer, rs2)
    length(rs2) == 0 && return empty!(rs1)

    rv2 = collect(rs2)
    sort!(rv2)

    intersect!(rs1, first(rv2):last(rv2))

    length(rv2) == 1 && return rs1

    if hasmethod(Iterators.reverse, Tuple{typeof(rv2)}) &&
       hasmethod(keys, Tuple{typeof(rv2)})
        rnext, rv2head = Iterators.peel(Iterators.reverse(rv2))
        for r in rv2head
            delete!(rs1, r+1:rnext-1)
            rnext = r
        end
    else
        rprev, rv2tail = Iterators.peel(rv2)
        for r in rv2tail
            delete!(rs1, rprev+1:r-1)
            rprev = r
        end
    end

    rs1
end

Base.intersect!(rs1::AbstractUnitRangesSortedContainer, rs2::AbstractSet{T}) where {T<:AbstractRange} = __intersect1!(rs1, rs2)
Base.intersect!(rs1::AbstractUnitRangesSortedContainer, rs2::AbstractVector{T}) where {T<:AbstractRange} = __intersect1!(rs1, rs2)
Base.intersect!(rs1::AbstractUnitRangesSortedContainer, rs2::Tuple{Vararg{AbstractRange}}) = __intersect1!(rs1, rs2)
function __intersect1!(rs1::AbstractUnitRangesSortedContainer, rs2)
    length(rs2) == 0 && return empty!(rs1)

    rv2 = collect(rs2)
    sort!(rv2)

    intersect!(rs1, first(rv2[firstindex(rv2)]):last(rv2[lastindex(rv2)]))

    length(rv2) == 1 && return rs1

    hasmethod(Iterators.reverse, Tuple{typeof(rv2)}) || throw(AssertionError("FIXME: TODO."))

    rnext, rv2head = Iterators.peel(Iterators.reverse(rv2))
    for r in rv2head
        delete!(rs1, last(r)+1:first(rnext)-1)
        rnext = r
    end


    rs1
end

function Base.intersect!(rs1::AbstractVector, rs2::AbstractUnitRangesSortedSet)
    rs = UnitRangesSortedSet(rs1)
    intersect!(rs, rs2)
    if eltype(rs1) <: AbstractRange
        resize!(rs1, length(rs))
        for (i, r) in enumerate(rs)
            rs1[i] = r
        end
    else
        resize!(rs1, sum(length(r) for r in rs))
        i = 0
        for r in rs, v in r
            rs1[i+=1] = v
        end
    end
    rs1
end
function Base.intersect!(rs1::AbstractSet, rs2::AbstractUnitRangesSortedSet)
    rs = UnitRangesSortedSet(rs1)
    intersect!(rs, rs2)
    empty!(rs1)
    if eltype(rs1) <: AbstractRange
        for r in rs
            push!(rs1, r)
        end
    else
        for r in rs, v in r
            push!(rs1, v)
        end
    end
    rs1
end

Base.intersect(rs1::AbstractUnitRangesSortedSet, rs2) = intersect!(copy(rs1), rs2)
Base.intersect(rs1::AbstractUnitRangesSortedSet, rs2::AbstractUnitRangesSortedSet) = intersect!(copy(rs1), rs2)
Base.intersect(rs1::AbstractSet, rs2::AbstractUnitRangesSortedSet) = _intersect2(rs1, rs2)  # AbstractSet ambiguity resolving
Base.intersect(rs1::AbstractVector, rs2::AbstractUnitRangesSortedSet) = _intersect2(rs1, rs2)
function _intersect2(rs1, rs2::AbstractUnitRangesSortedSet)
    rs = UnitRangesSortedSet(rs1)
    intersect!(rs, rs2)
    return convert(typeof(rs1), rs)
end


#
#  Aux functions
#

# https://github.com/JuliaLang/julia/issues/39952
basetype(::Type{T}) where T = Base.typename(T).wrapper
## https://discourse.julialang.org/t/how-do-a-i-get-a-type-stripped-of-parameters/73465/8
#basetype(::Type{T}) where T = eval(nameof(T))
## https://discourse.julialang.org/t/deparametrising-types/41939/4
#basetype(T::DataType) = T.name.wrapper
#basetype(T::UnionAll) = basetype(T.body)

datatypeshortname(x) = (T = typeof(x); string(basetype(T)) * "{" * string(eltype(eltype(T))) * "}")

Base.show(io::IO, ur::URSSIndexURange) = print(io, repr(first(ur)), ':', length(ur) == 0 ? "0" : "1" , ':', repr(last(ur)))

function Base.show(io::IO, ::MIME"text/plain", x::URSSIndexURange)
    len = length(x)
    print(io, len, "-element range ", first(x), ":", last(x), " in ", datatypeshortname(x.parent))
    if len != 0
        println(io, " containing:")
        show(IOContext(io, :typeinfo => eltype(x)), x)
    end
end

function Base.show(io::IOContext, x::URSSIndexURange)
    ranges = [getindex(x, s) for s in x]

    n = length(ranges)
    if isempty(ranges)
        return show(io, MIME("text/plain"), x)
    end
    limit = get(io, :limit, false)::Bool
    half_screen_rows = limit ? div(displaysize(io)[1] - 8, 2) : typemax(Int)
    pad = get_max_pad(x.parent)
    if !haskey(io, :compact)
        io = IOContext(io, :compact => true)
    end
    for k = eachindex(ranges)
        if k < half_screen_rows || k > length(ranges) - half_screen_rows
            print(io, "  ")
            if isassigned(ranges, Int(k))
                print(io, lpad(repr(first(ranges[k])), pad), ":", repr(last(ranges[k])))
            else
                print(io, Base.undef_ref_str)
            end
            k != length(ranges) && println(io)
        elseif k == half_screen_rows
            println(io, " "^(pad-1), "   \u22ee")
        end
    end
end

function Base.show(io::IO, x::T) where {T<:AbstractUnitRangesSortedSubSet}
    print(io, "subset(", datatypeshortname(x.parent), ", ", x.firstindex, ":", x.lastindex,  "): (")
    if length(x) > 0
        el1, restx = Iterators.peel(x)
        print(io, repr(first(el1)), ":", repr(last(el1)))
        for r in restx
            print(io, ", ", repr(first(r)), ":", repr(last(r)))
        end
    end
    print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", x::T) where {T<:AbstractUnitRangesSortedSubSet}
    len = length(x)
    print(io, len, "-element subset(", datatypeshortname(x.parent), ", ", x.firstindex, ":", x.lastindex, ")")
    if len != 0
        println(io, ":")
        show(IOContext(io, :typeinfo => eltype(x)), x)
    end
end


function Base.show(io::IO, x::T) where {T<:AbstractUnitRangesSortedSet}
    print(io, datatypeshortname(x), "(")
    if length(x) > 0
        el1, restx = Iterators.peel(x)
        print(io, repr(first(el1)), ":", repr(last(el1)))
        for r in restx
            print(io, ", ", repr(first(r)), ":", repr(last(r)))
        end
    end
    print(io, ")")
end

# derived from stdlib/SparseArrays/src/sparsevector.jl
function Base.show(io::IO, ::MIME"text/plain", x::T) where {T<:AbstractUnitRangesSortedSet}
    len = length(x)
    if len == 0
        println(io, datatypeshortname(x), "()")
    else
        println(io, datatypeshortname(x), " with ", len, len==1 ? " element:" : " elements:")
        show(IOContext(io, :typeinfo => eltype(x)), x)
    end
end

Base.show(io::IOContext, x::T) where {T<:AbstractUnitRangesSortedSubSet} = _display_AURSS(io, x)
Base.show(io::IOContext, x::T) where {T<:AbstractUnitRangesSortedSet} = _display_AURSS(io, x)

function _display_AURSS(io, x)
    ranges = [v for v in x]
    n = length(ranges)
    if isempty(ranges)
        return show(io, MIME("text/plain"), x)
    end
    limit = get(io, :limit, false)::Bool
    half_screen_rows = limit ? div(displaysize(io)[1] - 8, 2) : typemax(Int)
    pad = get_max_pad(x)
    if !haskey(io, :compact)
        io = IOContext(io, :compact => true)
    end
    for k = eachindex(ranges)
        if k < half_screen_rows || k > length(ranges) - half_screen_rows
            print(io, "  ")
            if isassigned(ranges, Int(k))
                print(io, lpad(repr(first(ranges[k])), pad), ":", repr(last(ranges[k])))
            else
                print(io, Base.undef_ref_str)
            end
            k != length(ranges) && println(io)
        elseif k == half_screen_rows
            println(io, " "^(pad-1), "   \u22ee")
        end
    end
end

function get_max_pad(rs::AbstractUnitRangesSortedSet)
    len = length(rs)
    pad = 0
    for (i,r) in enumerate(rs)
        if i < 100 || i > len - 100
            pad = max(pad, length(repr(first(r))), length(repr(last(r))))
        end
    end
    pad
end


end  # of module UnitRangesSortedSets
