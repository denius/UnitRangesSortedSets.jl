
module UnitRangesSortedSets
export AbstractUnitRangesSortedSet, AbstractSubUnitRangesSortedSet, AbstractUnitRangesSortedContainer
export UnitRangesSortedVector, UnitRangesSortedSet
export Sub0UnitRangesSortedSet, Sub1UnitRangesSortedSet, SubUnitRangesSortedSet, URSSIndexURange
export testfun_create, testfun_createSV, testfun_createVL, testfun_create_seq, testfun_create_dense, testfun_delete!,
       testfun_in, testfun_in_outer, testfun_in_rand, testfun_in_seq, testfun_nzgetindex, testfun_setindex!
export searchsortedrange, searchsortedfirstrange, searchsortedlastrange, getrange, getindex, beforestartindex, pastendindex
export subset


import Base: ForwardOrdering, Forward
const FOrd = ForwardOrdering

using DocStringExtensions
using DataStructures
import DataStructures: DataStructures.Tokens.IntSemiToken, DataStructures.SDMToken
using IterTools
using Setfield
using Random


# https://github.com/JuliaLang/julia/issues/39952
basetype(::Type{T}) where T = Base.typename(T).wrapper
## https://discourse.julialang.org/t/deparametrising-types/41939/4
#basetype(T::DataType) = T.name.wrapper
#basetype(T::UnionAll) = basetype(T.body)


@inline to_urange(::Type{TU}, l, r) where {TU<:UnitRange{Ti}} where Ti =
    TU(l, r)
@inline to_urange(::Type{TU}, II::AbstractRange) where {TU<:UnitRange{Ti}} where Ti =
    TU(first(II), last(II))
@inline to_urange(::Type{TU}, l, r) where {TU<:StepRange{Ti,Tst}} where {Ti,Tst} =
    TU(l, Tst(1), r)
@inline to_urange(::Type{TU}, II::AbstractRange) where {TU<:StepRange{Ti,Tst}} where {Ti,Tst} =
    TU(first(II), Tst(1), last(II))


@inline inferrangetype(::Type{Ti}) where {Ti<:Real} = UnitRange{Ti}
@inline function inferrangetype(::Type{Ti}) where Ti
    T = typeof(Ti(0):Ti(0))
    return T <: StepRange ? StepRange{Ti,Int8} : T
end


abstract type AbstractUnitRangesSortedSet{Ti,TU} <: AbstractSet{TU} end
abstract type AbstractUnitRangesSortedContainer{Ti,TU} <: AbstractUnitRangesSortedSet{Ti,TU} end
abstract type AbstractSubUnitRangesSortedSet{Ti,TU,P} <: AbstractUnitRangesSortedSet{Ti,TU} end

"""
$(TYPEDEF)
Immutable struct fields:
$(TYPEDFIELDS)
"""
struct Sub0UnitRangesSortedSet{Ti,TU,P,Tix} <: AbstractSubUnitRangesSortedSet{Ti,TU,P}
    #"Empty range `kstart:kstop` for simplicity."
    #singlerange::TU
    "The `<:AbstractUnitRangesSortedSet` which subset point to."
    parent::P
    "`kstart:kstop` is the bounds of subset."
    kstart::Ti
    kstop::Ti
    "Empty range `firstindex:lastindex` locates the point in `parent` where `kstart:kstop` should be inserted."
    firstindex::Tix
    lastindex::Tix
    "Pre `firstindex` position, used for acceleration of `iteration`."
    beforestartindex::Tix
    "Post `lastindex` position, used for acceleration of `iteration`."
    pastendindex::Tix
    "Length of subset."
    numranges::Int
end

"""
$(TYPEDEF)
Immutable struct fields:
$(TYPEDFIELDS)
"""
struct Sub1UnitRangesSortedSet{Ti,TU,P,Tix} <: AbstractSubUnitRangesSortedSet{Ti,TU,P}
    "The `<:AbstractUnitRangesSortedSet` which subset point to."
    parent::P
    "`kstart:kstop` is the bounds of subset."
    kstart::Ti
    kstop::Ti
    "Single continuous range `kstart:kstop` from `parent` stored for simplicity."
    singlerange::TU
    "Range `firstindex:lastindex` locates the point in `parent` where `kstart:kstop` should be inserted."
    firstindex::Tix
    lastindex::Tix
    "Pre `firstindex` position, used for acceleration of `iteration`."
    beforestartindex::Tix
    "Post `lastindex` position, used for acceleration of `iteration`."
    pastendindex::Tix
    "Length of subset."
    numranges::Int
end

"""
$(TYPEDEF)
Immutable struct fields:
$(TYPEDFIELDS)
"""
struct SubUnitRangesSortedSet{Ti,TU,P,Tix} <: AbstractSubUnitRangesSortedSet{Ti,TU,P}
    "The `<:AbstractUnitRangesSortedSet` which subset point to."
    parent::P
    "`kstart:kstop` is the bounds of subset."
    kstart::Ti
    kstop::Ti
    "Saved first range `parent[firstindex]` truncated to `kstart` from left."
    firstrange::TU
    "Saved last range `parent[lastindex]` truncated to `kstop` from right."
    lastrange::TU
    "Range `firstindex:lastindex` locates the indices diapason in `parent` where `kstart:kstop` should be inserted."
    firstindex::Tix
    lastindex::Tix
    "Pre `firstindex` position, used for acceleration of `iterate`."
    beforestartindex::Tix
    "Post `lastindex` position, used for acceleration of `iterate`."
    pastendindex::Tix
    "Length of subset."
    numranges::Int
end

subset(rs::P, II::AbstractRange) where {P<:AbstractSubUnitRangesSortedSet} = subset(rs.parent, II)
function subset(rs::P, II::AbstractRange) where {P<:AbstractUnitRangesSortedContainer{Ti,TU}} where {Ti,TU}
    I = to_urange(TU, II)
    ir = searchsortedrange(rs, I)
    if length(I) == 0 || length(ir) == 0
        if length(I) == 0 && length(ir) == 1
            ir = create_indexrange(rs, first(ir), regress(rs, first(ir)))
        end
        if length(I) != 0
            I = to_urange(TU, last(I), first(I))
        end
        return Sub0UnitRangesSortedSet{Ti,TU,P,typeof(first(ir))}(rs, first(I), last(I), first(ir), last(ir),
                                       regress(rs, first(ir)), advance(rs, last(ir)), length(ir))
    elseif length(ir) == 1
        singlerange = getindex(rs, first(ir))
        if first(singlerange) < first(I)
            singlerange = to_urange(TU, first(I), last(singlerange))
        end
        if last(I) < last(singlerange)
            singlerange = to_urange(TU, first(singlerange), last(I))
        end
        return Sub1UnitRangesSortedSet{Ti,TU,P,typeof(first(ir))}(rs, first(I), last(I), singlerange,
                                       first(ir), last(ir), regress(rs, first(ir)), advance(rs, last(ir)), length(ir))
    else
        firstrange = getindex(rs, first(ir))
        if first(firstrange) < first(I)
            firstrange = to_urange(TU, first(I), last(firstrange))
        end
        lastrange = getindex(rs, last(ir))
        if last(I) < last(lastrange)
            lastrange = to_urange(TU, first(lastrange), last(I))
        end
        return SubUnitRangesSortedSet{Ti,TU,P,typeof(first(ir))}(rs, first(I), last(I), firstrange, lastrange,
                                      first(ir), last(ir), regress(rs, first(ir)), advance(rs, last(ir)), length(ir))
    end
end

@inline Base.copy(rs::T) where {T<:AbstractSubUnitRangesSortedSet{Ti,TU}} where {Ti,TU} =
    intersect(rs.parent, to_urange(TU, rs.kstart, rs.kstop))

"""
Inserting zero, or negative length ranges does nothing.
$(TYPEDEF)
Mutable struct fields:
$(TYPEDFIELDS)
"""
mutable struct UnitRangesSortedVector{Ti,TU} <: AbstractUnitRangesSortedContainer{Ti,TU}
    "Index of last used range."
    lastusedrangeindex::Int
    "Storage for ranges starts,"
    rstarts::Vector{Ti}
    "and stops."
    rstops::Vector{Ti}
end

UnitRangesSortedVector{Ti,TU}() where {Ti,TU} = UnitRangesSortedVector{Ti,TU}(0, Vector{Ti}(undef, 0), Vector{Ti}(undef, 0))
UnitRangesSortedVector{Ti}() where {Ti} = UnitRangesSortedVector{Ti,inferrangetype(Ti)}()
function UnitRangesSortedVector(rs::AbstractUnitRangesSortedSet{Ti,TU}) where {Ti,TU}
    rstarts = Vector{Ti}(undef, length(rs))
    rstops = Vector{Ti}(undef, length(rs))
    for (i, r) in enumerate(rs)
        rstarts[i] = first(r)
        rstops[i] = last(r)
    end
    UnitRangesSortedVector{Ti,TU}(firstindex(rstarts) - 1, rstarts, rstops)
end


"""
$(TYPEDEF)
Mutable struct fields:
$(TYPEDFIELDS)
"""
mutable struct UnitRangesSortedSet{Ti,TU} <: AbstractUnitRangesSortedContainer{Ti,TU}
    "Index of last used range."
    lastusedrangeindex::IntSemiToken
    "Storage for ranges: the ket of Dict is the `first(range)`, and the value of Dict is the `last(range)`."
    ranges::SortedDict{Ti,Ti,FOrd}
end

function UnitRangesSortedSet{Ti,TU}() where {Ti,TU}
    ranges = SortedDict{Ti,Ti,FOrd}(Forward)
    UnitRangesSortedSet{Ti,TU}(beforestartsemitoken(ranges), ranges)
end
UnitRangesSortedSet{Ti}() where {Ti} = UnitRangesSortedSet{Ti,inferrangetype(Ti)}()
function UnitRangesSortedSet(rs::AbstractUnitRangesSortedSet{Ti,TU}) where {Ti,TU}
    ranges = SortedDict{Ti,Ti,FOrd}(Forward)
    for r in rs
        ranges[first(r)] = last(r)
    end
    UnitRangesSortedSet{Ti,TU}(beforestartsemitoken(ranges), ranges)
end


(::Type{T})(values::Union{AbstractVector, AbstractSet, Tuple}) where {T<:AbstractUnitRangesSortedContainer} =
    eltype(values) <: AbstractRange ? T{eltype(eltype(values))}(values) : T{eltype(values)}(values)

(::Type{T})(values::Union{AbstractVector, AbstractSet, Tuple}) where {T<:AbstractUnitRangesSortedContainer{Ti}} where {Ti} =
    T{inferrangetype(Ti)}(values)

function (::Type{T})(values::Union{AbstractVector, AbstractSet, Tuple}) where {T<:AbstractUnitRangesSortedContainer{Ti,TU}} where {Ti,TU}
    rs = T()
    for r in values
        push!(rs, r)
    end
    rs
end

Base.convert(::Type{<:UnitRangesSortedVector}, rs::UnitRangesSortedVector) = rs
Base.convert(::Type{<:UnitRangesSortedSet}, rs::UnitRangesSortedSet) = rs
function Base.convert(::Type{T}, rs::Union{UnitRangesSortedVector,UnitRangesSortedSet}) where {T<:AbstractVector{Tv}} where {Tv<:AbstractRange}
    V = T(undef, length(rs))
    for (i, r) in enumerate(rs)
        V[i] = to_urange(Tv, r)
    end
    V
end
function Base.convert(::Type{T}, rs::Union{UnitRangesSortedVector,UnitRangesSortedSet}) where {T<:AbstractVector{Tv}} where {Tv}
    V = T(undef, sum(length(r) for r in rs))
    i = 0
    for r in rs, v in r
        V[i+=1] = v
    end
    V
end
function Base.convert(::Type{T}, rs::Union{UnitRangesSortedVector,UnitRangesSortedSet}) where {T<:AbstractSet{Tv}} where {Tv<:AbstractRange}
    S = T()
    for r in rs
        push!(S, to_urange(Tv, r))
    end
    S
end
function Base.convert(::Type{T}, rs::Union{UnitRangesSortedVector,UnitRangesSortedSet}) where {T<:AbstractSet{Tv}} where {Tv}
    S = T()
    for r in rs, v in r
        push!(S, v)
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
@inline create_indexrange(parent::P, l, r) where {P<:UnitRangesSortedVector{Ti,TU}} where {Ti,TU} = UnitRange{Int}(l, r)
@inline create_indexrange(parent::P, l, r) where {P<:UnitRangesSortedSet{Ti,TU}} where {Ti,TU} = URSSIndexURange(parent, l, r)


@inline function URSSIndexURange(rs::P, l::Tix, r::Tix) where {P<:UnitRangesSortedVector, Tix}
    @boundscheck l < beforestartindex(rs) || l > pastendindex(rs) && return throw(BoundsError(rs, l))
    @boundscheck r < beforestartindex(rs) || r > pastendindex(rs) && return throw(BoundsError(rs, r))
    if l == pastendindex(rs) ||
       r == beforestartindex(rs) ||
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
    if l == pastendindex(rs) ||
       r == beforestartindex(rs) ||
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
                   {Pl<:UnitRangesSortedSet,Tixl,Pr<:UnitRangesSortedVector,Tixr} = false
@inline Base.:(==)(l::URSSIndexURange{Pl,Tixl}, r::URSSIndexURange{Pr,Tixr}) where
                   {Pl<:UnitRangesSortedVector,Tixl,Pr<:UnitRangesSortedSet,Tixr} = false
@inline Base.:(==)(l::URSSIndexURange{Pl,Tixl}, r::URSSIndexURange{Pr,Tixr}) where
                   {Pl<:UnitRangesSortedVector,Tixl,Pr<:UnitRangesSortedVector,Tixr} =
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

#@inline Base.collect(ur::URSSIndexURange{P,Tix}) where {P<:UnitRangesSortedVector,Tix} =
#    map(s->getindex(ur.parent, s), ur.start:ur.stop)
#@inline Base.collect(ur::URSSIndexURange{P,Tix}) where {P<:UnitRangesSortedSet,Tix} =
#    map(s->getindex(ur.parent, s), onlysemitokens(inclusive(ur.parent.ranges, ur.start, ur.stop)))

@inline DataStructures.advance(ur::URSSIndexURange, st::IntSemiToken) = advance((ur.parent, st))
@inline DataStructures.regress(ur::URSSIndexURange, st::IntSemiToken) = regress((ur.parent, st))
@inline beforestartindex(ur::URSSIndexURange) = beforestartsemitoken(ur.parent)
@inline pastendindex(ur::URSSIndexURange) = pastendsemitoken(ur.parent)

Base.length(rs::UnitRangesSortedVector) = length(rs.rstarts)
Base.length(rs::UnitRangesSortedSet) = length(rs.ranges)
Base.length(rs::AbstractSubUnitRangesSortedSet) = rs.numranges
Base.haslength(rs::AbstractUnitRangesSortedSet) = true
Base.hasfastin(rs::AbstractUnitRangesSortedSet) = true
Base.isempty(rs::AbstractUnitRangesSortedSet) = length(rs) == 0
Base.size(rs::AbstractUnitRangesSortedSet) = (length(rs),)
#Base.axes(rs::AbstractUnitRangesSortedSet) = (Base.OneTo(rs.n),)
Base.eltype(::AbstractUnitRangesSortedSet{Ti,TU}) where {Ti,TU} = TU
#Base.IndexStyle(::AbstractUnitRangesSortedSet) = IndexLinear()


function Base.collect(::Type{ElType}, rs::AbstractUnitRangesSortedSet{Ti,TU}) where {ElType,Ti,TU}
    T = basetype(TU){ElType}
    res = Vector{T}(undef, length(rs))
    i = 0
    for r in rs
        res[i+=1] = to_urange(T, ElType(first(r)), ElType(last(r)))
    end
    return res
end
function Base.collect(rs::AbstractUnitRangesSortedSet{Ti,TU}) where {Ti,TU}
    res = Vector{TU}(undef, length(rs))
    i = 0
    for r in rs
        res[i+=1] = to_urange(TU, r)
    end
    return res
end

@inline function indexcompare(rs::UnitRangesSortedVector, i, j)
    if i < j
        return -1
    elseif i > j
        return 1
    else
        return 0
    end
end
@inline function indexcompare(rs::AbstractSubUnitRangesSortedSet{Ti,TU,P}, i, j) where {Ti,TU,P<:UnitRangesSortedVector}
    if i < j
        return -1
    elseif i > j
        return 1
    else
        return 0
    end
end
@inline indexcompare(rs::UnitRangesSortedSet, i, j) = compare(rs.ranges, i, j)
@inline indexcompare(rs::AbstractSubUnitRangesSortedSet{Ti,TU,P}, i, j) where {Ti,TU,P<:UnitRangesSortedSet} =
    compare(rs.parent.ranges, i, j)


@inline function getrange(rs::AbstractUnitRangesSortedSet{Ti,TU}, i) where {Ti,TU}
    if in(i, rs)
        return to_urange(TU, first(i), last(i))
    else
        return nothing
    end
end

@inline getindex_tuple(rs::UnitRangesSortedVector, i) = (rs.rstarts[i], rs.rstops[i])
@inline getindex_tuple(rs::UnitRangesSortedSet, i::IntSemiToken) = tuple(deref((rs.ranges, i))...)
@inline getindex_tuple(rs::Sub0UnitRangesSortedSet, i) = (r = getindex(rs, i); tuple(first(r), last(r)))
@inline getindex_tuple(rs::Sub1UnitRangesSortedSet, i) = (r = getindex(rs, i); tuple(first(r), last(r)))
@inline getindex_tuple(rs::SubUnitRangesSortedSet, i) = (r = getindex(rs, i); tuple(first(r), last(r)))
@inline Base.getindex(rs::UnitRangesSortedVector{Ti,TU}, i) where {Ti,TU} =
    to_urange(TU, rs.rstarts[i], rs.rstops[i])
@inline Base.getindex(rs::UnitRangesSortedSet{Ti,TU}, i::IntSemiToken) where {Ti,TU} =
    to_urange(TU, deref((rs.ranges, i))...)
@inline Base.getindex(rs::UnitRangesSortedSet{Ti,TU}, t::SDMToken) where {Ti,TU} =
    to_urange(TU, deref(t)...)
@inline function Base.getindex(ur::URSSIndexURange, i)
    @boundscheck indexcompare(ur.parent, first(ur), i) != 1 &&
                 indexcompare(ur.parent, i, last(ur)) != 1 || throw(BoundsError(ur, i))
    getindex(ur.parent, i)
end
@inline Base.getindex(rs::Sub0UnitRangesSortedSet{Ti,TU}, i) where {Ti,TU} =
    (throw(BoundsError(rs)); to_urange(TU, typemin(Ti), typemin(Ti)))
@inline function Base.getindex(rs::Sub1UnitRangesSortedSet, i)
    @boundscheck i == rs.firstindex || throw(BoundsError(rs, i))
    rs.singlerange
end
@inline function Base.getindex(rs::SubUnitRangesSortedSet{Ti,TU,P}, i) where {Ti,TU,P<:UnitRangesSortedVector}
    @boundscheck beforestartindex(rs) < i < pastendindex(rs) || throw(BoundsError(rs, i))
    if i == rs.firstindex
        return rs.firstrange
    elseif i < rs.lastindex
        return getindex(rs.parent, i)
    else
        return rs.lastrange
    end
end
@inline function Base.getindex(rs::SubUnitRangesSortedSet, i)
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

@inline getindex_rangestart(rs::UnitRangesSortedVector, i) = rs.rstarts[i]
@inline getindex_rangestop(rs::UnitRangesSortedVector, i) = rs.rstops[i]
@inline getindex_rangestart(rs::UnitRangesSortedSet, i) = deref_key((rs.ranges, i))
@inline getindex_rangestop(rs::UnitRangesSortedSet, i) = deref_value((rs.ranges, i))
@inline getindex_rangestart(rs::AbstractSubUnitRangesSortedSet, i) = first(getindex(rs, i))
@inline getindex_rangestop(rs::AbstractSubUnitRangesSortedSet, i) = last(getindex(rs, i))

@inline Base.firstindex(rs::UnitRangesSortedVector) = firstindex(rs.rstarts)
@inline Base.firstindex(rs::UnitRangesSortedSet) = startof(rs.ranges)
@inline Base.firstindex(rs::AbstractSubUnitRangesSortedSet) = rs.firstindex
@inline Base.lastindex(rs::UnitRangesSortedVector) = lastindex(rs.rstarts)
@inline Base.lastindex(rs::UnitRangesSortedSet) = lastindex(rs.ranges)
@inline Base.lastindex(rs::AbstractSubUnitRangesSortedSet) = rs.lastindex
@inline Base.first(rs::UnitRangesSortedVector{Ti,TU}) where {Ti,TU} =
    to_urange(TU, rs.rstarts[1], rs.rstops[1])
@inline Base.first(rs::UnitRangesSortedSet{Ti,TU}) where {Ti,TU} =
    to_urange(TU, deref((rs.ranges, startof(rs.ranges)))...)
@inline Base.first(rs::Sub0UnitRangesSortedSet) = getindex(rs, firstindex(rs))
@inline Base.first(rs::Sub1UnitRangesSortedSet) = rs.singlerange
@inline Base.first(rs::SubUnitRangesSortedSet) = getindex(rs, firstindex(rs))
@inline Base.last(rs::UnitRangesSortedVector{Ti,TU}) where {Ti,TU} =
    to_urange(TU, rs.rstarts[end], rs.rstops[end])
@inline Base.last(rs::UnitRangesSortedSet{Ti,TU}) where {Ti,TU} =
    to_urange(TU, deref((rs.ranges, lastindex(rs.ranges)))...)
@inline Base.last(rs::Sub0UnitRangesSortedSet) = getindex(rs, lastindex(rs))
@inline Base.last(rs::Sub1UnitRangesSortedSet) = rs.singlerange
@inline Base.last(rs::SubUnitRangesSortedSet) = getindex(rs, lastindex(rs))

@inline beforestartindex(rs::UnitRangesSortedVector) = firstindex(rs.rstarts) - 1
@inline beforestartindex(rs::UnitRangesSortedSet) = beforestartsemitoken(rs.ranges)
@inline beforestartindex(rs::AbstractSubUnitRangesSortedSet) = rs.beforestartindex
@inline pastendindex(rs::UnitRangesSortedVector) = lastindex(rs.rstarts) + 1
@inline pastendindex(rs::UnitRangesSortedSet) = pastendsemitoken(rs.ranges)
@inline pastendindex(rs::AbstractSubUnitRangesSortedSet) = rs.pastendindex

@inline DataStructures.advance(rs::UnitRangesSortedVector, state) = state + 1
@inline DataStructures.advance(rs::UnitRangesSortedSet, state) = advance((rs.ranges, state))
@inline DataStructures.advance(rs::AbstractSubUnitRangesSortedSet, state) = advance(rs.parent, state)
@inline DataStructures.regress(rs::UnitRangesSortedVector, state) = state - 1
@inline DataStructures.regress(rs::UnitRangesSortedSet, state) = regress((rs.ranges, state))
@inline DataStructures.regress(rs::AbstractSubUnitRangesSortedSet, state) = regress(rs.parent, state)

# Derived from julia/base/sort.jl, for increasing sorted vectors.
function bisectionsearchlast(V::AbstractVector, val)
    lo = 0
    hi = length(V) + 1
    @inbounds while lo < hi - 1
        m = lo + ((hi - lo) >>> 0x01)
        if val < V[m]
            hi = m
        else
            lo = m
        end
    end
    return lo
end
function bisectionsearchlast(V::AbstractVector, val, lo::T, hi::T) where T
    u = T(1)
    lo = lo - u
    hi = hi + u
    @inbounds while lo < hi - u
        m = lo + ((hi - lo) >>> 0x01)
        if val < V[m]
            hi = m
        else
            lo = m
        end
    end
    return lo
end

function bisectionsearchfirst(V::AbstractVector, val)
    lo = 0
    hi = length(V) + 1
    @inbounds while lo < hi - 1
        m = lo + ((hi - lo) >>> 0x01)
        if V[m] < val
            lo = m
        else
            hi = m
        end
    end
    return hi
end
function bisectionsearchfirst(V::AbstractVector, val, lo::T, hi::T) where T
    u = T(1)
    lo = lo - u
    hi = hi + u
    @inbounds while lo < hi - u
        m = lo + ((hi - lo) >>> 0x01)
        if V[m] < val
            lo = m
        else
            hi = m
        end
    end
    return hi
end


"Returns index of range in which, or after, `k` is placed."
@inline searchsortedlastrange(rs::UnitRangesSortedVector{Ti}, k) where {Ti} = bisectionsearchlast(rs.rstarts, Ti(k))
#@inline searchsortedlastrange(rs::UnitRangesSortedVector{Ti}, k) where {Ti} = searchsortedlast(rs.rstarts, Ti(k); lt=<)
@inline searchsortedlastrange(rs::UnitRangesSortedSet{Ti}, k) where {Ti} = searchsortedlast(rs.ranges, Ti(k))
@inline searchsortedlastrange(rs::Sub0UnitRangesSortedSet, k) = rs.beforestartindex
@inline searchsortedlastrange(rs::Sub1UnitRangesSortedSet{Ti}, k) where {Ti} =
    Ti(k) >= rs.kstart ? rs.firstindex : rs.beforestartindex
@inline searchsortedlastrange(rs::SubUnitRangesSortedSet{Ti,TU,P,Tix}, k) where {Ti,TU,P<:UnitRangesSortedVector,Tix} =
    bisectionsearchlast(rs.parent.rstarts, Ti(k), rs.firstindex, rs.lastindex)
@inline function searchsortedlastrange(rs::SubUnitRangesSortedSet{Ti,TU,P,Tix}, k) where {Ti,TU,P<:UnitRangesSortedSet,Tix}
    st = searchsortedlast(rs.parent.ranges, Ti(k))
    if compare(rs.parent.ranges, st, rs.beforestartindex) == 1 && compare(rs.parent.ranges, st, rs.pastendindex) == -1
        return st
    else
        return beforestartindex(rs)
    end
end

"Returns index of range in which, or before, `k` is placed."
@inline searchsortedfirstrange(rs::UnitRangesSortedVector{Ti}, k) where {Ti} = bisectionsearchfirst(rs.rstops, Ti(k))
@inline function searchsortedfirstrange(rs::UnitRangesSortedSet{Ti}, k) where {Ti}
    st = searchsortedlastrange(rs, Ti(k))
    if st != beforestartindex(rs) && in(Ti(k), getindex(rs, st))
        return st
    else
        return advance(rs, st)
    end
end
@inline searchsortedfirstrange(rs::Sub0UnitRangesSortedSet, k) = rs.pastendindex
@inline searchsortedfirstrange(rs::Sub1UnitRangesSortedSet{Ti}, k) where {Ti} = Ti(k) <= rs.kstop ? rs.lastindex : rs.pastendindex
@inline searchsortedfirstrange(rs::SubUnitRangesSortedSet{Ti,TU,P,Tix}, k) where {Ti,TU,P<:UnitRangesSortedVector,Tix} =
    bisectionsearchfirst(rs.parent.rstops, Ti(k), rs.firstindex, rs.lastindex)
@inline function searchsortedfirstrange(rs::SubUnitRangesSortedSet{Ti,TU,P,Tix}, k) where {Ti,TU,P<:UnitRangesSortedSet,Tix}
    st = searchsortedlastrange(rs, Ti(k))
    if st != beforestartindex(rs) && in(Ti(k), getindex(rs, st))
        return st
    else
        return advance(rs, st)
    end
end

"Returns range of `rs` indexes which coincide or concluded in `I` range."
@inline searchsortedrange(rs::AbstractUnitRangesSortedContainer, I::AbstractRange) =
    create_indexrange(rs, searchsortedfirstrange(rs, first(I)), searchsortedlastrange(rs, last(I)))
@inline searchsortedrange(rs::AbstractSubUnitRangesSortedSet, I::AbstractRange) =
    create_indexrange(rs.parent, searchsortedfirstrange(rs, first(I)), searchsortedlastrange(rs, last(I)))

"Returns indexes of range in `rs` in which `k` may be inserted. Or negative range in the case of `k` is
 between `rs` ranges, and indices of resulted range is the indexes of that neighbors."
@inline function searchsortedrange(rs::AbstractUnitRangesSortedContainer{Ti}, k) where {Ti}
    st = searchsortedlastrange(rs, Ti(k))
    if st != beforestartindex(rs) && in(Ti(k), getindex(rs, st))
        return create_indexrange(rs, st, st)
    else
        return create_indexrange(rs, advance(rs, st), st)
    end
end


@inline index_status(rs::UnitRangesSortedVector, st) = 0 <= st <= length(rs)+1 ? 1 : 0
@inline index_status(rs::UnitRangesSortedSet, st) = status((rs.ranges, st))

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

@inline Base.iterate(rs::Sub0UnitRangesSortedSet{Ti,TU}, state = 1) where {Ti,TU} =
    # for type inference
    rs.kstart > rs.kstop ? nothing : (to_urange(TU, typemin(Ti), typemin(Ti)), state)
@inline Base.iterate(rrs::Base.Iterators.Reverse{T}, state = 1) where {T<:Sub0UnitRangesSortedSet{Ti,TU}} where {Ti,TU} =
    rrs.itr.kstart > rrs.itr.kstop ? nothing : (to_urange(TU, typemin(Ti), typemin(Ti)), state)


@inline Base.iterate(rs::Sub1UnitRangesSortedSet) = (rs.singlerange, rs.firstindex)
@inline Base.iterate(rs::Sub1UnitRangesSortedSet, state) = nothing
@inline Base.iterate(rrs::Base.Iterators.Reverse{T}) where {T<:Sub1UnitRangesSortedSet} =
    (rrs.itr.singlerange, rrs.itr.lastindex)
@inline Base.iterate(rrs::Base.Iterators.Reverse{T}, state) where {T<:Sub1UnitRangesSortedSet} = nothing


@inline function Base.iterate(rs::AbstractSubUnitRangesSortedSet{Ti,TU,P},
                              state = (rs.firstrange, rs.firstindex)) where {Ti,TU,P<:UnitRangesSortedVector}
    val, state2 = state
    if state2 < rs.lastindex
        state2 += 1
        return (val, (getindex(rs.parent, state2), state2))
    elseif state2 == rs.lastindex
        state2 += 1
        return (rs.lastrange, (rs.lastrange, state2))
    else
        return nothing
    end
end
@inline function Base.iterate(rrs::Base.Iterators.Reverse{T}, state = (rrs.itr.lastrange, rrs.itr.lastindex)) where
                              {T<:AbstractSubUnitRangesSortedSet{Ti,TU,P}} where {Ti,TU,P<:UnitRangesSortedVector}
    val, state2 = state
    if rrs.itr.firstindex < state2
        state2 -= 1
        return (val, (getindex(rrs.itr.parent, state2), state2))
    elseif rrs.itr.firstindex == state2
        state2 -= 1
        return (rrs.itr.firstrange, (rrs.itr.firstrange, state2))
    else
        return nothing
    end
end


@inline function Base.iterate(rs::AbstractSubUnitRangesSortedSet{Ti,TU,P},
                              state = (rs.firstrange, rs.firstindex)) where {Ti,TU,P<:UnitRangesSortedSet}
    val, state2 = state
    if compare(rs.parent.ranges, state2, rs.lastindex) == -1
        state2 = advance((rs.parent.ranges, state2))
        return (val, (getindex(rs.parent, state2), state2))
    elseif compare(rs.parent.ranges, state2, rs.lastindex) == 0
        state2 = advance((rs.parent.ranges, state2))
        return (rs.lastrange, (rs.lastrange, state2))
    else
        return nothing
    end
end
@inline function Base.iterate(rrs::Base.Iterators.Reverse{T}, state = (rrs.itr.lastrange, rrs.itr.lastindex)) where
                              {T<:AbstractSubUnitRangesSortedSet{Ti,TU,P}} where {Ti,TU,P<:UnitRangesSortedSet}
    val, state2 = state
    if compare(rrs.itr.parent.ranges, rrs.itr.firstindex, state2) == -1
        state2 = regress((rrs.itr.parent.ranges, state2))
        return (val, (getindex(rrs.itr.parent, state2), state2))
    elseif compare(rrs.itr.parent.ranges, rrs.itr.firstindex, state2) == 0
        state2 = regress((rrs.itr.parent.ranges, state2))
        return (rrs.itr.firstrange, (rrs.itr.firstrange, state2))
    else
        return nothing
    end
end



@inline function Base.iterate(rs::AbstractUnitRangesSortedSet, state = firstindex(rs))
    if state != pastendindex(rs)
        return (getindex(rs, state), advance(rs, state))
    else
        return nothing
    end
end
@inline function Base.iterate(rrs::Base.Iterators.Reverse{T}, state = lastindex(rrs.itr)) where {T<:AbstractUnitRangesSortedSet}
    if state != beforestartindex(rrs.itr)
        return (getindex(rrs.itr, state), regress(rrs.itr, state))
    else
        return nothing
    end
end

#@inline function Base.iterate(rs::UnitRangesSortedSet, state = SDMToken((rs.ranges, startof(rs.ranges))))
#    if state[2] != pastendsemitoken(state[1])
#        return (getindex(rs, state), SDMToken((state[1], advance(state))))
#    else
#        return nothing
#    end
#end

@inline function Base.iterate(ur::URSSIndexURange, state = (first(ur), 0))
    st, i = state #[1], state[2]
    if i < length(ur)
        return (st, (advance(ur.parent, st), i + 1))
    else
        return nothing
    end
end


@inline function eachindexiterate(rs::AbstractUnitRangesSortedSet, state = firstindex(rs))
    if state != pastendindex(rs)
        return (state, advance(rs, state))
    else
        return nothing
    end
end
@inline function eachindexiterate(rrs::Base.Iterators.Reverse{T}, state = lastindex(rrs.itr)) where {T<:AbstractUnitRangesSortedSet}
    if state != beforestartindex(rrs.itr)
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
#Base.IteratorEltype(::Type{EachIndexIterator{It}}) where {It} = Base.EltypeUnknown()
Base.IteratorSize(::Type{<:EachIndexIterator}) = Base.HasShape{1}()
Base.length(it::EachIndexIterator) = length(it.itr)
Base.size(it::EachIndexIterator) = (length(it.itr),)
Iterators.reverse(it::EachIndexIterator) = EachIndexIterator(Iterators.reverse(it.itr))


Base.eachindex(rs::UnitRangesSortedVector) = eachindex(rs.rstarts)
Base.eachindex(rs::UnitRangesSortedSet) = eachindexiterator(rs)
Base.eachindex(rs::T) where {T<:AbstractSubUnitRangesSortedSet} = eachindexiterator(rs)

#
# Assignments
#


@inline Base.in(II::AbstractRange, rs::AbstractUnitRangesSortedSet{Ti,TU}) where {Ti,TU} = _in(to_urange(TU, II), rs)
@inline Base.in(II::TU, rs::AbstractUnitRangesSortedSet{Ti,TU}) where {Ti,TU} = _in(II, rs)

@inline function _in(I::TU, rs::AbstractUnitRangesSortedSet{Ti,TU}) where {Ti,TU}
    # fast check for cached range index
    if (st = rs.lastusedrangeindex) != beforestartindex(rs)
        if issubset(I, getindex(rs, st))
            return true
        end
    end
    # cached range' index miss (or index is not stored), thus try search
    Iposition = searchsortedrange(rs, I)
    rs.lastusedrangeindex = last(Iposition)
    if length(Iposition) != 1
        return false
    elseif issubset(I, getindex(rs, first(Iposition)))
        return true
    else
        return false
    end
end

@inline function Base.in(key, rs::AbstractUnitRangesSortedSet{Ti}) where Ti
    k = Ti(key)
    # fast check for cached range index
    if (st = rs.lastusedrangeindex) != beforestartindex(rs)
        r_start, r_stop = getindex_tuple(rs, st)
        if r_start <= k <= r_stop
            return true
        end
    end
    # cached range index miss (or index not stored), thus try search
    st = searchsortedlastrange(rs, k)
    if st != beforestartindex(rs)  # `k` is not before the start of first range
        r_start, r_stop = getindex_tuple(rs, st)
        if k <= r_stop  # is `k` inside of range
            rs.lastusedrangeindex = st
            return true
        end
    end
    rs.lastusedrangeindex = beforestartindex(rs)
    return false
end

@inline Base.in(II::AbstractRange, su::Sub0UnitRangesSortedSet) = false
@inline Base.in(key, su::Sub0UnitRangesSortedSet) = false
@inline Base.in(II::AbstractRange, su::Sub1UnitRangesSortedSet{Ti,TU}) where {Ti,TU} = issubset(to_urange(TU, II), su.singlerange)
@inline Base.in(II::TU, su::Sub1UnitRangesSortedSet{Ti,TU}) where {Ti,TU} = issubset(II, su.singlerange)
@inline Base.in(key, su::Sub1UnitRangesSortedSet{Ti}) where {Ti} = in(Ti(key), su.singlerange)

@inline Base.in(II::AbstractRange, su::SubUnitRangesSortedSet{Ti,TU}) where {Ti,TU} = _in(to_urange(TU, II), su)
@inline Base.in(II::TU, su::SubUnitRangesSortedSet{Ti,TU}) where {Ti,TU} = _in(II, su)
@inline function _in(I::TU, su::SubUnitRangesSortedSet{Ti,TU}) where {Ti,TU}
    # fast check for cached range index
    if (st = su.parent.lastusedrangeindex) != beforestartindex(su.parent)
        if issubset(I, getindex(su, st))
            return true
        end
    end
    # cached range' index miss (or index is not stored), thus try search
    Iposition = searchsortedrange(su, I)
    su.parent.lastusedrangeindex = last(Iposition)
    if length(Iposition) != 1
        return false
    elseif issubset(I, getindex(su, first(Iposition)))
        return true
    else
        return false
    end
end

@inline function Base.in(key, su::SubUnitRangesSortedSet{Ti}) where {Ti}
    k = Ti(key)
    # fast check for cached range index
    if (st = su.parent.lastusedrangeindex) != beforestartindex(su.parent)
        r_start, r_stop = getindex_tuple(su, st)
        if r_start <= k <= r_stop
            return true
        end
    end
    # cached range index miss (or index not stored), thus try search
    st = searchsortedlastrange(su, k)
    if st != beforestartindex(su)  # `k` is not before the start of first range
        r_start, r_stop = getindex_tuple(su, st)
        if r_start <= k <= r_stop  # is `k` inside of range
            su.parent.lastusedrangeindex = st
            return true
        end
    end
    su.parent.lastusedrangeindex = beforestartindex(su.parent)
    return false
end


#@inline Base.haskey(rs::AbstractUnitRangesSortedSet, key) = in(key, rs)

function Base.push!(rs::UnitRangesSortedVector{Ti,TU}, II::Union{AbstractVector,AbstractSet,NTuple}) where {Ti,TU}
    for r in II
        push!(rs, r)
    end
    rs
end


@inline Base.push!(rs::UnitRangesSortedVector{Ti,TU}, II::AbstractRange) where {Ti,TU} = _push!(rs, to_urange(TU, II))
@inline Base.push!(rs::UnitRangesSortedVector{Ti,TU}, II::TU) where {Ti,TU} = _push!(rs, II)

function _push!(rs::UnitRangesSortedVector{Ti,TU}, I::TU) where {Ti,TU}

    if length(I) < 2
        for i in I
            push!(rs, i)
        end
        return rs
    end

    rs.lastusedrangeindex = beforestartindex(rs)

    Iposition = searchsortedrange(rs, I)

    # `I` already exist in `rs` and the same as the corresponding range in `rs`, nothing to do
    if length(Iposition) == 1 && I == getindex(rs, first(Iposition))
        return rs

    # delete all inside `I` range in `rs`
    else
        delete!(rs, I)
    end

    # get neighbors pointers
    # ATTENTION: Note: the range will be zero-length, thus
    # `first(Iposition)` will be point to range in `rs` on right side to `I` and
    # `last(Iposition)` will be point to range in `rs` on left side to `I`.
    Iposition = searchsortedrange(rs, first(I))
    @boundscheck @assert Iposition == searchsortedrange(rs, last(I)) "FIXME: Something went wrong."

    Iposition_left = last(Iposition)
    Iposition_right = first(Iposition)

    # `I` is adjoined in both sides with `rs` ranges, thus join them all
    if Iposition_left != beforestartindex(rs) && getindex_rangestop(rs, Iposition_left) + 1 == first(I) &&
       Iposition_right != pastendindex(rs) && getindex_rangestart(rs, Iposition_right) - 1 == last(I)
        rs.rstops[Iposition_left] = getindex_rangestop(rs, Iposition_right)
        deleteat!(rs.rstarts, Iposition_right)
        deleteat!(rs.rstops, Iposition_right)

    # `I` is adjoin with `rs` range on left side, thus append to left range
    elseif Iposition_left != beforestartindex(rs) && getindex_rangestop(rs, Iposition_left) + 1 == first(I) &&
           Iposition_right != pastendindex(rs) && getindex_rangestart(rs, Iposition_right) - 1 != last(I)
        rs.rstops[Iposition_left] = last(I)

    # `I` is adjoin with `rs` range on right side, thus prepend to right range
    elseif Iposition_left != beforestartindex(rs) && getindex_rangestop(rs, Iposition_left) + 1 != first(I) &&
           Iposition_right != pastendindex(rs) && getindex_rangestart(rs, Iposition_right) - 1 == last(I)
        rs.rstarts[Iposition_right] = first(I)

    # `I` is separate from both sides, insert it
    elseif Iposition_left != beforestartindex(rs) && getindex_rangestop(rs, Iposition_left) + 1 != first(I) &&
           Iposition_right != pastendindex(rs) && getindex_rangestart(rs, Iposition_right) - 1 != last(I)
        insert!(rs.rstarts, Iposition_right, first(I))
        insert!(rs.rstops, Iposition_right, last(I))

    # `I` is first range in `rs` and adjoin with range on right side, thus prepend `I` to right range
    elseif Iposition_left == beforestartindex(rs) &&
           Iposition_right != pastendindex(rs) && getindex_rangestart(rs, Iposition_right) - 1 == last(I)
        rs.rstarts[Iposition_right] = first(I)

    # `I` is separate first range in `rs`, insert it
    elseif Iposition_left == beforestartindex(rs) &&
           Iposition_right != pastendindex(rs) && getindex_rangestart(rs, Iposition_right) - 1 != last(I)
        insert!(rs.rstarts, Iposition_right, first(I))
        insert!(rs.rstops, Iposition_right, last(I))

    # `I` is last range in `rs` and adjoin with range on left side, thus append `I` to left range
    elseif Iposition_left != beforestartindex(rs) && getindex_rangestop(rs, Iposition_left) + 1 == first(I) &&
           Iposition_right == pastendindex(rs)
        rs.rstops[Iposition_left] = last(I)

    # `I` is separate last range in `rs`, insert it
    elseif Iposition_left != beforestartindex(rs) && getindex_rangestop(rs, Iposition_left) + 1 != first(I) &&
           Iposition_right == pastendindex(rs)
        insert!(rs.rstarts, Iposition_right, first(I))
        insert!(rs.rstops, Iposition_right, last(I))

    # `rs` is empty
    elseif Iposition_left == beforestartindex(rs) && Iposition_right == pastendindex(rs)
        insert!(rs.rstarts, 1, first(I))
        insert!(rs.rstops, 1, last(I))

    else
        throw(AssertionError("FIXME: Something went wrong."))
    end

    return rs
end

function Base.push!(rs::UnitRangesSortedSet{Ti,TU}, II::Union{AbstractVector,AbstractSet,NTuple}) where {Ti,TU}
    for r in II
        push!(rs, r)
    end
    rs
end

@inline Base.push!(rs::UnitRangesSortedSet{Ti,TU}, II::AbstractRange) where {Ti,TU} = _push!(rs, to_urange(TU, II))
@inline Base.push!(rs::UnitRangesSortedSet{Ti,TU}, II::TU) where {Ti,TU} = _push!(rs, II)
function _push!(rs::UnitRangesSortedSet{Ti,TU}, I::TU) where {Ti,TU}

    if length(I) < 2
        for i in I
            push!(rs, i)
        end
        return rs
    end

    rs.lastusedrangeindex = beforestartindex(rs)

    Iposition = searchsortedrange(rs, I)

    # `I` already exist in `rs` and the same as the corresponding range in `rs`, nothing to do
    if length(Iposition) == 1 && I == getindex(rs, first(Iposition))
        return rs

    # delete all inside `I` range in `rs`
    else
        delete!(rs, I)
    end

    # get neighbors pointers
    # ATTENTION: Note: the range will be zero-length, thus
    # `first(Iposition)` will be point to range in `rs` on right side to `I` and
    # `last(Iposition)` will be point to range in `rs` on left side to `I`.
    Iposition = searchsortedrange(rs, first(I))
    @boundscheck @assert Iposition == searchsortedrange(rs, last(I)) "FIXME: Something went wrong."

    Iposition_left = last(Iposition)
    Iposition_right = first(Iposition)
    Iposition_left_rangestop = Iposition_left != beforestartindex(rs) ? getindex_rangestop(rs, Iposition_left) :
                                                                         typemin(Ti)
    Iposition_right_rangestart = Iposition_right != pastendindex(rs) ? getindex_rangestart(rs, Iposition_right) :
                                                                        typemax(Ti)

    # `I` is adjoined in both sides with `rs` ranges, thus join them all
    if Iposition_left != beforestartindex(rs) && Iposition_left_rangestop + 1 == first(I) &&
       Iposition_right != pastendindex(rs) && Iposition_right_rangestart - 1 == last(I)
        rs.ranges[Iposition_left] = getindex_rangestop(rs, Iposition_right)
        delete!((rs.ranges, Iposition_right))

    # `I` is adjoin with `rs` range on left side, thus append to left range
    elseif Iposition_left != beforestartindex(rs) && Iposition_left_rangestop + 1 == first(I) &&
           Iposition_right != pastendindex(rs) && Iposition_right_rangestart - 1 != last(I)
        rs.ranges[Iposition_left] = last(I)

    # `I` is adjoin with `rs` range on right side, thus prepend to right range
    elseif Iposition_left != beforestartindex(rs) && Iposition_left_rangestop + 1 != first(I) &&
           Iposition_right != pastendindex(rs) && Iposition_right_rangestart - 1 == last(I)
        rs.ranges[first(I)] = getindex_rangestop(rs, Iposition_right)
        delete!((rs.ranges, Iposition_right))

    # `I` is separate from both sides, insert it
    elseif Iposition_left != beforestartindex(rs) && Iposition_left_rangestop + 1 != first(I) &&
           Iposition_right != pastendindex(rs) && Iposition_right_rangestart - 1 != last(I)
        rs.ranges[first(I)] = last(I)

    # `I` is first range in `rs` and adjoin with range on right side, thus prepend `I` to right range
    elseif Iposition_left == beforestartindex(rs) &&
           Iposition_right != pastendindex(rs) && Iposition_right_rangestart - 1 == last(I)
        rs.ranges[first(I)] = getindex_rangestop(rs, Iposition_right)
        delete!((rs.ranges, Iposition_right))

    # `I` is separate first range in `rs`, insert it
    elseif Iposition_left == beforestartindex(rs) &&
           Iposition_right != pastendindex(rs) && Iposition_right_rangestart - 1 != last(I)
        rs.ranges[first(I)] = last(I)

    # `I` is last range in `rs` and adjoin with range on left side, thus append `I` to left range
    elseif Iposition_left != beforestartindex(rs) && Iposition_left_rangestop + 1 == first(I) &&
           Iposition_right == pastendindex(rs)
        rs.ranges[Iposition_left] = last(I)

    # `I` is separate last range in `rs`, insert it
    elseif Iposition_left != beforestartindex(rs) && Iposition_left_rangestop + 1 != first(I) &&
           Iposition_right == pastendindex(rs)
        rs.ranges[first(I)] = last(I)

    # `rs` is empty
    elseif Iposition_left == beforestartindex(rs) && Iposition_right == pastendindex(rs)
        rs.ranges[first(I)] = last(I)

    else
        throw(AssertionError("FIXME: Something went wrong."))
    end

    return rs
end

@inline function check_exist_and_update(rs, k)

    # fast check for cached range index
    if (st = rs.lastusedrangeindex) != beforestartindex(rs)
        r_start, r_stop = getindex_tuple(rs, st)
        if r_start <= k <= r_stop
            return nothing
        end
    end

    st = searchsortedlastrange(rs, k)
    #
    #sstatus = status((rs.ranges, st))
    @boundscheck if index_status(rs, st) == 0 # invalid range index
        throw(BoundsError(rs, k))
    end

    # check the index exist and update its data
    if st != beforestartindex(rs)  # `k` is not before the first index
        if k <= getindex_rangestop(rs, st)
            rs.lastusedrangeindex = st
            return nothing
        end
    end

    return st
end


function Base.push!(rs::UnitRangesSortedVector{Ti,TU}, key) where {Ti,TU}
    k = Ti(key)

    st = check_exist_and_update(rs, k)
    st == nothing && return rs

    if length(rs) == 0
        push!(rs.rstarts, k)
        push!(rs.rstops, k)
        rs.lastusedrangeindex = 1
        return rs
    end

    if st == beforestartindex(rs)  # `k` is before the first range
        if rs.rstarts[1] - k > 1  # there is will be gap in indices after inserting
            pushfirst!(rs.rstarts, k)
            pushfirst!(rs.rstops, k)
        else  # prepend to first range
            rs.rstarts[1] = k
        end
        rs.lastusedrangeindex = 1
        return rs
    end

    r = getindex(rs, st)

    if k >= rs.rstarts[end]  # `k` is after the last range start
        if k > last(r) + 1  # there is will be the gap in indices after inserting
            push!(rs.rstarts, k)
            push!(rs.rstops, k)
        else  # just append to last range
            rs.rstops[st] += 1
        end
        rs.lastusedrangeindex = lastindex(rs)
        return rs
    end

    # `k` is somewhere between indices
    stnext = advance(rs, st)
    rnext = getindex(rs, stnext)

    if first(rnext) - last(r) == 2  # join ranges
        rs.rstarts[st] = first(r)
        rs.rstops[st] = last(rnext)
        deleteat!(rs.rstarts, stnext)
        deleteat!(rs.rstops, stnext)
        rs.lastusedrangeindex = st
    elseif k - last(r) == 1  # append to left range
        rs.rstops[st] += 1
        rs.lastusedrangeindex = st
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


function Base.push!(rs::UnitRangesSortedSet{Ti,TU}, key) where {Ti,TU}
    k = Ti(key)

    st = check_exist_and_update(rs, k)
    st == nothing && return rs

    if length(rs) == 0
        rs.ranges[k] = k
        rs.lastusedrangeindex = firstindex(rs)
        return rs
    end

    if st == beforestartindex(rs)  # `k` is before the first index
        stnext = advance(rs, st)
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

    r = getindex(rs, st)

    if k > last(last(rs)) # `k` is after the last range end index
        if last(r) + 1 < k  # there is will be the gap in indices after inserting
            rs.ranges[k] = k
        else  # just append to last range
            rs.ranges[st] = last(r)+1
        end
        rs.lastusedrangeindex = lastindex(rs)
        return rs
    end

    rs.lastusedrangeindex = beforestartindex(rs)

    # `k` is somewhere between indices
    stnext = advance(rs, st)
    rnext = getindex(rs, stnext)

    if first(rnext) - last(r) == 2  # join ranges
        rs.ranges[st] = last(rnext)
        delete!((rs.ranges, stnext))
    elseif k - last(r) == 1  # append to left range
        rs.ranges[st] = k
        rs.lastusedrangeindex = st
    elseif first(rnext) - k == 1  # prepend to right range
        rs.ranges[k] = last(rnext)
        delete!((rs.ranges, stnext))
    else  # insert single element range
        rs.ranges[k] = k
    end

    return rs

end


function Base.delete!(rs::UnitRangesSortedVector{Ti,TU}, II::T) where {Ti,TU, T<:Union{AbstractVector,AbstractSet,NTuple}}
    if hasmethod(Iterators.reverse, Tuple{T})
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

@inline Base.delete!(rs::UnitRangesSortedVector{Ti,TU}, II::AbstractRange) where {Ti,TU} = _delete!(rs, to_urange(TU, II))
@inline Base.delete!(rs::UnitRangesSortedVector{Ti,TU}, II::TU) where {Ti,TU} = _delete!(rs, II)
function _delete!(rs::UnitRangesSortedVector{Ti,TU}, I::TU) where {Ti,TU}
    length(I) == 0 && return rs

    Iposition = searchsortedrange(rs, I)
    length(Iposition) == 0 && return rs

    rs.lastusedrangeindex = beforestartindex(rs)

    # delete all concluded ranges
    todelete1 = getindex_rangestart(rs, first(Iposition)) < first(I) ? first(Iposition) + 1 : first(Iposition)
    todelete2 = getindex_rangestop(rs, last(Iposition)) > last(I) ? last(Iposition) - 1 : last(Iposition)
    splice!(rs.rstarts, todelete1:todelete2)
    splice!(rs.rstops, todelete1:todelete2)

    # remains the side ranges to delete
    # Note: there is not possible `beforestartindex` or `pastendindex` in `Iposition`
    Iposition = searchsortedrange(rs, I)
    length(Iposition) == 0 && return rs

    Iposition_rangestart = getindex_rangestart(rs, first(Iposition))
    Iposition_rangestop = getindex_rangestop(rs, last(Iposition))

    # `I` intersects only one range in `rs`
    if length(Iposition) == 1

        # inside one range, thus split it
        if Iposition_rangestart < first(I) && last(I) < Iposition_rangestop
            insert!(rs.rstarts, advance(rs, first(Iposition)), last(I) + 1)
            insert!(rs.rstops, first(Iposition), first(I) - 1)

        # `I` intersects with or inside in one range from left side, thus shrink from left side
        elseif Iposition_rangestart >= first(I) && last(I) < Iposition_rangestop
            rs.rstarts[first(Iposition)] = last(I) + 1

        # `I` intersects with or inside in one range from right side, thus shrink from right side
        # inside one range, ended to left side, thus shrink from right side
        elseif Iposition_rangestart < first(I) && last(I) >= Iposition_rangestop
            rs.rstops[first(Iposition)] = first(I) - 1

        else
            throw(AssertionError("FIXME: Something went wrong."))
        end

    # All remaining cases are with two ranges from `rs`,
    # and both side ranges only partly intersects with `I`.
    elseif length(Iposition) == 2

        # shrink second (right) range from the left side
        rs.rstarts[last(Iposition)] = last(I) + 1

        # shrink first (left) range from the right side
        rs.rstops[first(Iposition)] = first(I) - 1

    else
        throw(AssertionError("FIXME: Something went wrong."))
    end

    return rs
end


function Base.delete!(rs::UnitRangesSortedSet{Ti,TU}, II::Union{AbstractVector,AbstractSet,NTuple})  where {Ti,TU}
    for r in II
        delete!(rs, r)
    end
    rs
end

@inline Base.delete!(rs::UnitRangesSortedSet{Ti,TU}, II::AbstractRange) where {Ti,TU} = _delete!(rs, to_urange(TU, II))
@inline Base.delete!(rs::UnitRangesSortedSet{Ti,TU}, II::TU) where {Ti,TU} = _delete!(rs, II)
function _delete!(rs::UnitRangesSortedSet{Ti,TU}, I::TU) where {Ti,TU}
    length(I) == 0 && return rs

    Iposition = searchsortedrange(rs, I)
    length(Iposition) == 0 && return rs

    rs.lastusedrangeindex = beforestartindex(rs)

    # delete all concluded ranges
    todelete1 = getindex_rangestart(rs, first(Iposition)) < first(I) ? advance(rs, first(Iposition)) : first(Iposition)
    todelete2 = getindex_rangestop(rs, last(Iposition)) > last(I) ? regress(rs, last(Iposition)) : last(Iposition)
    for st in onlysemitokens(inclusive(rs.ranges, todelete1, todelete2))
        delete!((rs.ranges, st))
    end

    # remains the side ranges to delete
    # Note: there is not possible `beforestartindex` or `pastendindex` in `Iposition`
    Iposition = searchsortedrange(rs, I)
    length(Iposition) == 0 && return rs

    Iposition_rangestart = getindex_rangestart(rs, first(Iposition))
    Iposition_rangestop = getindex_rangestop(rs, last(Iposition))

    # `I` intersects only one range in `rs`
    if length(Iposition) == 1

        # inside one range, thus split it
        if Iposition_rangestart < first(I) && last(I) < Iposition_rangestop
            rs.ranges[last(I) + 1] = Iposition_rangestop
            rs.ranges[first(Iposition)] = first(I) - 1

        # `I` intersects with or inside in one range from left side, thus shrink from left side
        elseif Iposition_rangestart >= first(I) && last(I) < Iposition_rangestop
            rs.ranges[last(I) + 1] = Iposition_rangestop
            delete!((rs.ranges, first(Iposition)))

        # `I` intersects with or inside in one range from right side, thus shrink from right side
        # inside one range, ended to left side, thus shrink from right side
        elseif Iposition_rangestart < first(I) && last(I) >= Iposition_rangestop
            rs.ranges[first(Iposition)] = first(I) - 1

        else
            throw(AssertionError("FIXME: Something went wrong."))
        end

    # All remaining cases are with two ranges from `rs`,
    # and both side ranges only partly intersects with `I`.
    elseif length(Iposition) == 2

        # shrink second (right) range from the left side
        rs.ranges[last(I) + 1] = Iposition_rangestop
        delete!((rs.ranges, last(Iposition)))


        # shrink first (left) range from the right side
        rs.ranges[first(Iposition)] = first(I) - 1

    else
        throw(AssertionError("FIXME: Something went wrong."))
    end

    return rs
end


function Base.delete!(rs::UnitRangesSortedVector{Ti,TU}, key) where {Ti,TU}
    k = Ti(key)

    length(rs) == 0 && return rs

    st = searchsortedlastrange(rs, k)

    if st == beforestartindex(rs)  # `k` is before first range
        return rs
    end

    r = getindex(rs, st)

    if k > last(r)  # `k` is outside of range
        return rs
    end

    if last(r) - first(r) + 1 == 1
        deleteat!(rs.rstarts, st)
        deleteat!(rs.rstops, st)
    elseif k == last(r)  # last index in range
        rs.rstops[st] -= 1
    elseif k == first(r)  # first index in range
        rs.rstarts[st] += 1
    else
        insert!(rs.rstarts, advance(rs, st), k+1)
        insert!(rs.rstops, advance(rs, st), last(r))
        rs.rstops[st] = k-1
    end

    rs.lastusedrangeindex = beforestartindex(rs)

    return rs
end


function Base.delete!(rs::UnitRangesSortedSet{Ti,TU}, key) where {Ti,TU}
    k = Ti(key)

    length(rs) == 0 && return rs

    st = searchsortedlastrange(rs, k)

    if st == beforestartindex(rs)  # `k` is before first index
        return rs
    end

    r = getindex(rs, st)

    if k > last(r)  # `k` is outside of range
        return rs
    end

    if last(r) - first(r) + 1 == 1
        delete!((rs.ranges, st))
    elseif k == last(r)  # last index in range
        rs.ranges[st] = last(r) - 1
    elseif k == first(r)  # first index in range
        delete!((rs.ranges, st))
        rs.ranges[k+1] = last(r)
    else
        rs.ranges[st] = k - 1
        rs.ranges[k+1] = last(r)
    end

    rs.lastusedrangeindex = beforestartindex(rs)

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


# TODO: setdiff, symdiff, eachindex (SomeEachindexIterator), map
# setdiff already exist in julia/base/abstractset.jl


Base.empty(rs::T) where {T<:AbstractUnitRangesSortedContainer} = T()

function Base.empty!(rs::UnitRangesSortedVector)
    empty!(rs.rstarts)
    empty!(rs.rstops)
    rs.lastusedrangeindex = beforestartindex(rs)
    rs
end
function Base.empty!(rs::AbstractUnitRangesSortedContainer)
    empty!(rs.ranges)
    rs.lastusedrangeindex = beforestartindex(rs)
    rs
end

Base.copy(rs::T) where {T<:UnitRangesSortedVector} = T(rs.lastusedrangeindex, copy(rs.rstarts), copy(rs.rstops))
Base.copy(rs::T) where {T<:UnitRangesSortedSet} = T(rs.lastusedrangeindex, packcopy(rs.ranges))

#Base.deepcopy(rs::AbstractUnitRangesSortedContainer) = copy(rs)



Base.union(rs::AbstractUnitRangesSortedSet, rss...) = union!(copy(rs), rss...)
Base.union(rs::AbstractUnitRangesSortedSet, rs2) = union!(copy(rs), rs2)

@inline Base.union!(rs::AbstractUnitRangesSortedContainer, rss...) = union!(union!(rs, rss[1]), Base.tail(rss)...)
function Base.union!(rs::AbstractUnitRangesSortedContainer, rs2)
    issubset(rs2, rs) && return rs
    for r in rs2
        push!(rs, r)
    end
    return rs
end

@inline Base.union!(rs::AbstractUnitRangesSortedContainer, r::AbstractRange) = push!(rs, r)


Base.setdiff(rs::AbstractUnitRangesSortedSet, rss...) = setdiff!(copy(rs), rss...)
Base.setdiff(rs::AbstractUnitRangesSortedSet, rs2) = setdiff!(copy(rs), rs2)

@inline Base.setdiff!(rs::AbstractUnitRangesSortedContainer, rss...) = setdiff!(setdiff!(rs, rss[1]), Base.tail(rss)...)
function Base.setdiff!(rs::AbstractUnitRangesSortedContainer, rs2)
    if hasmethod(Iterators.reverse, Tuple{typeof(rs2)})
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

Base.symdiff(rs::AbstractUnitRangesSortedSet, sets...) = symdiff!(empty(rs), rs, sets...)
Base.symdiff(rs::AbstractUnitRangesSortedSet, s) = symdiff!(copy(rs), s)

#function Base.symdiff!(rs::AbstractUnitRangesSortedSet, itrs...)
#    for x in itrs
#        symdiff!(rs, x)
#    end
#    return rs
#end

@inline Base.symdiff!(rs::AbstractUnitRangesSortedContainer, rss...) = symdiff!(symdiff!(rs, rss[1]), Base.tail(rss)...)
function Base.symdiff!(rs1::AbstractUnitRangesSortedContainer, rs2)
    for r in rs2
        ss = subset(rs1, r)
        vv = collect(ss)
        if length(vv) == 0
            push!(rs1, r)
        elseif length(vv) == 1
            delete!(rs1, vv[1])
            push!(rs1, first(r):first(vv[1])-1)
            push!(rs1, last(vv[1])+1:last(r))
        else
            push!(rs1, first(r):first(vv[1])-1)
            push!(rs1, last(vv[end])+1:last(r))
            rnext, rhead = Iterators.peel(Iterators.reverse(vv))
            delete!(rs1, rnext)
            for rr in rhead
                delete!(rs1, rr)
                push!(rs1, last(rr)+1:first(rnext)-1)
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
function Base.issubset(rs1::AbstractUnitRangesSortedSet, rs2::Union{AbstractSet,AbstractVector,AbstractRange,Tuple})
    for r1 in rs1
        issubset(r1, rs2) || return false
    end
    return true
end
function Base.issubset(rs1::AbstractUnitRangesSortedSet, rs2::Union{AbstractSet{T},AbstractVector{T},NTuple{N,T}}) where {T<:AbstractRange,N}
    for r1 in rs1
        findfirst(s->issubset(r1, s), rs2) != nothing || return false
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

@inline Base.intersect!(rs::AbstractUnitRangesSortedContainer{Ti,TU}, II::AbstractRange) where {Ti,TU} =
    __intersect!(rs, to_urange(TU,II))
@inline Base.intersect!(rs::AbstractUnitRangesSortedContainer{Ti,TU}, II::TU) where {Ti,TU} = __intersect!(rs, II)
function __intersect!(rs::AbstractUnitRangesSortedContainer{Ti,TU}, I::TU) where {Ti,TU}
    length(I) == 0 && return empty!(rs)
    delete!(rs, last(I) + 1:getindex_rangestop(rs, lastindex(rs)))
    delete!(rs, getindex_rangestart(rs, firstindex(rs)):first(I) - 1)
    rs
end

function Base.intersect!(rs1::AbstractUnitRangesSortedContainer, rs2::AbstractUnitRangesSortedSet)
    length(rs2) == 0 && return empty!(rs1)

    # cut both sides up to `rs2` indices
    intersect!(rs1, getindex_rangestart(rs2, firstindex(rs2)):getindex_rangestop(rs2, lastindex(rs2)))

    length(rs2) == 1 && return rs1

    if hasmethod(Iterators.reverse, Tuple{typeof(rs2)})
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
Base.intersect!(rs1::AbstractUnitRangesSortedContainer, rs2::AbstractSet) = _intersect!(rs1, rs2)
Base.intersect!(rs1::AbstractUnitRangesSortedContainer, rs2::AbstractVector) = _intersect!(rs1, rs2)
function _intersect!(rs1::AbstractUnitRangesSortedContainer, rs2)
    length(rs2) == 0 && return empty!(rs1)

    rv2 = collect(rs2)
    sort!(rv2)

    intersect!(rs1, first(rv2):last(rv2))

    length(rv2) == 1 && return rs1

    if hasmethod(Iterators.reverse, Tuple{typeof(rv2)})
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

function Base.intersect!(rs1::AbstractUnitRangesSortedContainer, rs2::Union{AbstractSet{T},AbstractVector{T},NTuple{N,T}}) where {T<:AbstractRange,N}
    length(rs2) == 0 && return empty!(rs1)

    rv2 = collect(rs2)
    sort!(rv2)

    intersect!(rs1, first(rv2[firstindex(rv2)]):last(rv2[lastindex(rv2)]))

    length(rv2) == 1 && return rs1

    if hasmethod(Iterators.reverse, Tuple{typeof(rv2)})
        rnext, rv2head = Iterators.peel(Iterators.reverse(rv2))
        for r in rv2head
            delete!(rs1, last(r)+1:first(rnext)-1)
            rnext = r
        end
    else
        rprev, rv2tail = Iterators.peel(rv2)
        for r in rv2tail
            delete!(rs1, last(rprev)+1:first(r)-1)
            rprev = r
        end
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
# AbstractSet ambiguity resolving
Base.intersect(rs1::AbstractSet, rs2::AbstractUnitRangesSortedSet) = _intersect(rs1, rs2)
Base.intersect(rs1::AbstractVector, rs2::AbstractUnitRangesSortedSet) = _intersect(rs1, rs2)
function _intersect(rs1, rs2::AbstractUnitRangesSortedSet)
    rs = UnitRangesSortedSet(rs1)
    intersect!(rs, rs2)
    return convert(typeof(rs1), rs)
end


#
#  Aux functions
#

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

function Base.show(io::IO, x::T) where {T<:AbstractSubUnitRangesSortedSet}
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

function Base.show(io::IO, ::MIME"text/plain", x::T) where {T<:AbstractSubUnitRangesSortedSet}
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
    print(io, datatypeshortname(x), "()")
    if len != 0
        println(io, ":")
        show(IOContext(io, :typeinfo => eltype(x)), x)
    end
end

Base.show(io::IOContext, x::T) where {T<:AbstractSubUnitRangesSortedSet} = _display_AURSS(io, x)
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

#
#  Testing functions
#

function testfun_create(T::Type, n = 500_000, density = 0.9)
    rs = T()
    Random.seed!(1234)
    randseq = randsubseq(1:n, density)
    rrandseq = shuffle(randseq)
    for (j,i) in enumerate(rrandseq)
        push!(rs, i)
    end
    for i in randseq
        in(i, rs) || println("Not coincide on index $i")
    end
    sum(length(r) for r in rs) == length(randseq) || println("Lost some indices")
    rs
end
function testfun_createSV(T::Type, n = 500_000, m = 5, density = 0.9)
    rs = T(m,n)
    Random.seed!(1234)
    for i in shuffle(randsubseq(1:n, density))
        for j = 1:m
            rs[i,j] = rand()
        end
    end
    rs
end
function testfun_createVL(T::Type, n = 500_000, density = 0.9)
    rs = T(n)
    Random.seed!(1234)
    for i in shuffle(randsubseq(1:n, density))
        rs[i] = rand(rand(0:7))
    end
    rs
end

function testfun_create_seq(T::Type, n = 500_000, density = 0.9)
    rs = T()
    Random.seed!(1234)
    randseq = randsubseq(1:n, density)
    for i in randseq
        push!(rs, i)
    end
    for i in randseq
        in(i, rs) || println("Not coincide on index $i")
    end
    sum(length(r) for r in rs) == length(randseq) || println("Lost some indices")
    rs
end

function testfun_create_dense(T::Type, n = 500_000, nchunks = 800, density = 0.95)
    rs = T()
    chunklen = max(1, floor(Int, n / nchunks))
    Random.seed!(1234)
    for i = 0:nchunks-1
        len = floor(Int, chunklen*density + randn() * chunklen * min(0.1, (1.0-density), density))
        len = max(1, min(chunklen-2, len))
        for j = 1:len
            push!(rs, i*chunklen + j)
        end
    end
    rs
end


function testfun_delete!(rs)
    Random.seed!(1234)
    indices = shuffle([x for r in rs for x in r])
    for i in indices
        delete!(rs, i)
    end
    length(rs) == 0 || println("Not all elements deleted")
    rs
end


function testfun_in(rs)
    start = first(first(rs))
    stop = last(last(rs))
    Random.seed!(1234)
    indices = shuffle(collect(start:stop))
    I = 0
    for i in indices
        in(i, rs) && (I += 1)
    end
    I
end

function testfun_in_outer(rs, idx)
    I = 0
    for i in idx
        in(i, rs) && (I += 1)
    end
    I
end
function testfun_in_rand(rs)
    len = sum(length(r) for r in rs)
    Random.seed!(1234)
    I = 0
    for j = 1:len
        i = rand(1:len)
        in(i, rs) && (I += 1)
    end
    I
end

function testfun_in_seq(rs)
    start = first(first(rs))
    stop = last(last(rs))
    I = 0
    for i in start:stop
        in(i, rs) && (I += 1)
    end
    I
end

function testfun_nzgetindex(sv)
    S = 0.0
    for i in nzindices(sv)
        S += sv[i]
    end
    (0, S)
end

function testfun_setindex!(sv)
    for i in nzindices(sv)
        sv[i] = 0.0
    end
end




end  # of module UnitRangesSortedSets
