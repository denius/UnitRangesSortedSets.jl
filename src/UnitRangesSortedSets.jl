
module UnitRangesSortedSets
export AbstractUnitRangesSortedSet, UnitRangesSortedVector, UnitRangesSortedSet, SubUnitRangesSortedSet, URSSUnitRange
export testfun_create, testfun_createSV, testfun_createVL, testfun_create_seq, testfun_create_dense, testfun_delete!,
       testfun_in, testfun_in_outer, testfun_in_rand, testfun_in_seq, testfun_nzgetindex, testfun_setindex!
export searchsortedrange, searchsortedfirstrange, searchsortedlastrange, findrange, getindex, beforestartindex, pastendindex


import Base: ForwardOrdering, Forward
const FOrd = ForwardOrdering

using DocStringExtensions
using DataStructures
using IterTools
using Setfield
using Random


abstract type AbstractUnitRangesSortedSet{Ti} <: AbstractSet{Ti} end

struct Sub0UnitRangesSortedSet{Ti,P,Tidx} <: AbstractUnitRangesSortedSet{Ti}
    data::UnitRange{Ti}
    parent::P
    start::Ti
    stop::Ti
    firstindex::Tidx
    lastindex::Tidx
    firstindexnext::Tidx
    lastindexprev::Tidx
    beforestartindex::Tidx
    pastendindex::Tidx
    numranges::Int
end

struct Sub1UnitRangesSortedSet{Ti,P,Tidx} <: AbstractUnitRangesSortedSet{Ti}
    data::UnitRange{Ti}
    parent::P
    start::Ti
    stop::Ti
    firstindex::Tidx
    lastindex::Tidx
    firstindexnext::Tidx
    lastindexprev::Tidx
    beforestartindex::Tidx
    pastendindex::Tidx
    numranges::Int
end

struct SubUnitRangesSortedSet{Ti,P,Tidx} <: AbstractUnitRangesSortedSet{Ti}
    parent::P
    start::Ti
    stop::Ti
    firstindex::Tidx
    lastindex::Tidx
    firstindexnext::Tidx
    lastindexprev::Tidx
    beforestartindex::Tidx
    pastendindex::Tidx
    numranges::Int
end

function Base.intersect(rs::Tp, II::AbstractRange) where {Tp<:AbstractUnitRangesSortedSet{Ti}} where Ti
    I = UnitRange{Ti}(first(II), last(II))
    ir = searchsortedrange(rs, I)
    if length(I) == 0 && length(ir) == 1
        ir = UnitRange(first(ir) + 1, last(ir))
    end
    next = advance(rs, first(ir))
    prev = regress(rs, last(ir))
    beforestart = first(ir) != beforestartindex(rs) ? regress(rs, first(ir)) : first(ir)
    pastend = last(ir) != pastendindex(rs) ? advance(rs, last(ir)) : last(ir)
    if length(I) == 0 || length(ir) == 0
        rdata = UnitRange(first(I), last(I))
        return Sub0UnitRangesSortedSet{Ti,Tp,typeof(first(ir))}(rdata, rs, first(I), last(I), first(ir), last(ir),
                                                               next, prev, beforestart, pastend, length(ir))
    elseif length(ir) == 1
        rdata = getindex(rs, first(ir))
        if first(rdata) < first(I)
            rdata = UnitRange(first(I), last(rdata))
        end
        if last(I) < last(rdata)
            rdata = UnitRange(first(rdata), last(I))
        end
        return Sub1UnitRangesSortedSet{Ti,Tp,typeof(first(ir))}(rdata, rs, first(I), last(I),
                                                               first(ir), last(ir), next, prev, beforestart, pastend, length(ir))
    else
        return SubUnitRangesSortedSet{Ti,Tp,typeof(first(ir))}(rs, first(I), last(I), first(ir), last(ir),
                                                              next, prev, beforestart, pastend, length(ir))
    end
end

Base.copy(rs::T) where {T<:Union{Sub0UnitRangesSortedSet,Sub1UnitRangesSortedSet,SubUnitRangesSortedSet}} =
    T(rs.parent, rs.start, rs.stop, rs.firstindex, rs.lastindex, rs.firstindexnext, rs.lastindexprev,
      rs.beforestartindex, rs.pastendindex, rs.numranges)
Base.deepcopy(rs::T) where {T<:Union{Sub0UnitRangesSortedSet,Sub1UnitRangesSortedSet,SubUnitRangesSortedSet}} =
    T(deepcopy(rs.parent), rs.start, rs.stop, rs.firstindex, rs.lastindex, rs.firstindexnext, rs.lastindexprev,
      rs.beforestartindex, rs.pastendindex, rs.numranges)

"""
Inserting zero, or negative length ranges does nothing.
$(TYPEDEF)
Mutable struct fields:
$(TYPEDFIELDS)
"""
mutable struct UnitRangesSortedVector{Ti} <: AbstractUnitRangesSortedSet{Ti}
    "Index of last used range"
    lastusedrangeindex::Int
    "Storage for ranges"
    rstarts::Vector{Ti}
    rstops::Vector{Ti}
end

UnitRangesSortedVector{Ti}() where {Ti} = UnitRangesSortedVector{Ti}(0, Vector{Ti}(undef, 0), Vector{Ti}(undef, 0))
function UnitRangesSortedVector(rs::AbstractUnitRangesSortedSet{Ti}) where {Ti}
    rstarts = Vector{Ti}(undef, length(rs))
    rstops = Vector{Ti}(undef, length(rs))
    for (i, r) in enumerate(rs)
        rstarts[i] = first(r)
        rstops[i] = last(r)
    end
    UnitRangesSortedVector{Ti}(firstindex(rstarts) - 1, rstarts, rstops)
end


"""
$(TYPEDEF)
Mutable struct fields:
$(TYPEDFIELDS)
"""
mutable struct UnitRangesSortedSet{Ti} <: AbstractUnitRangesSortedSet{Ti}
    "Index of last used range"
    lastusedrangeindex::DataStructures.Tokens.IntSemiToken
    "Storage for ranges: the ket of Dict is `first(range)`, and the value of Dict is the `last(range)`"
    ranges::SortedDict{Ti,Ti,FOrd}
end

function UnitRangesSortedSet{Ti}() where {Ti}
    ranges = SortedDict{Ti,Ti,FOrd}(Forward)
    UnitRangesSortedSet{Ti}(beforestartsemitoken(ranges), ranges)
end
function UnitRangesSortedSet(rs::AbstractUnitRangesSortedSet{Ti}) where {Ti}
    ranges = SortedDict{Ti,Ti,FOrd}(Forward)
    for r in rs
        ranges[first(r)] = last(r)
    end
    UnitRangesSortedSet{Ti}(beforestartsemitoken(ranges), ranges)
end


function (::Type{T})(values::Union{AbstractVector, AbstractSet, Tuple}) where {T<:AbstractUnitRangesSortedSet}
    if eltype(values) <: AbstractRange
        rs = T{eltype(eltype(values))}()
    else
        rs = T{eltype(values)}()
    end
    for r in values
        push!(rs, r)
    end
    rs
end
function (::Type{T})(values::Union{AbstractVector, AbstractSet, Tuple}) where {T<:AbstractUnitRangesSortedSet{Ti}} where Ti
    rs = T()
    for r in values
        push!(rs, r)
    end
    rs
end


"Type of `UnitRange` for `DataStructures.Tokens.IntSemiToken`."
struct URSSUnitRange{Trs,Tidx} <: AbstractUnitRange{Tidx}
    start::Tidx
    stop::Tidx
    parent::Trs
    len::Int
end


@inline function URSSUnitRange(rs::Trs, l::Tidx, r::Tidx) where {Trs<:UnitRangesSortedVector, Tidx}
    l < beforestartindex(rs) || l > pastendindex(rs) && return throw(BoundsError(rs, l))
    r < beforestartindex(rs) || r > pastendindex(rs) && return throw(BoundsError(rs, r))
    if l == pastendindex(rs) ||
       r == beforestartindex(rs) ||
       l > r
        len = 0
    else
        len = r - l + 1
    end
    URSSUnitRange{Trs,Tidx}(l, r, rs, len)
end

@inline function URSSUnitRange(rs::Trs, l::Tidx, r::Tidx) where {Trs<:UnitRangesSortedSet, Tidx<:DataStructures.Tokens.IntSemiToken}
    status((rs.ranges, l)) == 0 && return throw(KeyError(l))
    status((rs.ranges, r)) == 0 && return throw(KeyError(r))
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
    URSSUnitRange{Trs,Tidx}(l, r, rs, len)
end

#@inline Base.length(ur::URSSUnitRange{Trs,Tidx,Tl}) where {Trs<:UnitRangesSortedVector, Tidx, Tl} =
#    ur.len == 0 ? 0 : ur.stop - ur.start + 1
#@inline Base.length(ur::URSSUnitRange{Trs,Tidx,Tl}) where {Trs<:UnitRangesSortedSet, Tidx, Tl} = ur.len
@inline Base.length(ur::URSSUnitRange) = ur.len
@inline Base.step(rs::URSSUnitRange) = 1

@inline Base.:(==)(l::URSSUnitRange{Trsl,Tidxl}, r::URSSUnitRange{Trsr,Tidxr}) where
                   {Trsl<:UnitRangesSortedSet,Tidxl,Trsr<:UnitRangesSortedVector,Tidxr} = false
@inline Base.:(==)(l::URSSUnitRange{Trsl,Tidxl}, r::URSSUnitRange{Trsr,Tidxr}) where
                   {Trsl<:UnitRangesSortedVector,Tidxl,Trsr<:UnitRangesSortedSet,Tidxr} = false
@inline Base.:(==)(l::URSSUnitRange{Trsl,Tidxl}, r::URSSUnitRange{Trsr,Tidxr}) where
                   {Trsl<:UnitRangesSortedVector,Tidxl,Trsr<:UnitRangesSortedVector,Tidxr} =
    l.parent === r.parent &&
    l.start == r.start &&
    l.stop == r.stop
@inline Base.:(==)(l::URSSUnitRange{Trsl,Tidxl}, r::URSSUnitRange{Trsr,Tidxr}) where
                   {Trsl<:UnitRangesSortedSet,Tidxl,Trsr<:UnitRangesSortedSet,Tidxr} =
    l.parent === r.parent &&
    compare(l.parent.ranges, l.start, r.start) == 0 &&
    compare(l.parent.ranges, l.stop, r.stop) == 0

@inline Base.firstindex(ur::URSSUnitRange) = ur.start
@inline Base.lastindex(ur::URSSUnitRange) = ur.stop
@inline DataStructures.advance(ur::URSSUnitRange, st::DataStructures.Tokens.IntSemiToken) = advance((ur.parent, st))
@inline DataStructures.regress(ur::URSSUnitRange, st::DataStructures.Tokens.IntSemiToken) = regress((ur.parent, st))
@inline beforestartindex(ur::URSSUnitRange) = beforestartsemitoken(ur.parent)
@inline pastendindex(ur::URSSUnitRange) = pastendsemitoken(ur.parent)

@inline function Base.iterate(ur::URSSUnitRange, state = (firstindex(ur), 1))
    st, i = state[1], state[2]
    if i <= length(ur)
        return (st, (advance(ur.parent, st), i + 1))
    else
        return nothing
    end
end

#@inline Base.collect(ur::URSSUnitRange{Trs,Tidx}) where {Trs<:UnitRangesSortedVector,Tidx} =
#    map(s->getindex(ur.parent, s), ur.start:ur.stop)
#@inline Base.collect(ur::URSSUnitRange{Trs,Tidx}) where {Trs<:UnitRangesSortedSet,Tidx} =
#    map(s->getindex(ur.parent, s), onlysemitokens(inclusive(ur.parent.ranges, ur.start, ur.stop)))

Base.length(rs::UnitRangesSortedVector) = length(rs.rstarts)
Base.length(rs::UnitRangesSortedSet) = length(rs.ranges)
Base.length(rs::Union{Sub0UnitRangesSortedSet,Sub1UnitRangesSortedSet,SubUnitRangesSortedSet}) = rs.numranges
Base.haslength(rs::AbstractUnitRangesSortedSet) = true
Base.hasfastin(rs::AbstractUnitRangesSortedSet) = true
Base.isempty(rs::AbstractUnitRangesSortedSet) = length(rs) == 0
Base.size(rs::AbstractUnitRangesSortedSet) = (length(rs),)
#Base.axes(rs::AbstractUnitRangesSortedSet) = (Base.OneTo(rs.n),)
Base.eltype(::AbstractUnitRangesSortedSet{Ti}) where {Ti} = Ti
#Base.IndexStyle(::AbstractUnitRangesSortedSet) = IndexLinear()


function Base.collect(::Type{ElType}, rs::AbstractUnitRangesSortedSet) where ElType
    res = Vector{UnitRange{ElType}}(undef, length(rs))
    i = 0
    for r in rs
        res[i+=1] = UnitRange(ElType(first(r)), ElType(last(r)))
    end
    return res
end
function Base.collect(rs::AbstractUnitRangesSortedSet{Ti}) where Ti
    res = Vector{UnitRange{Ti}}(undef, length(rs))
    i = 0
    for r in rs
        res[i+=1] = UnitRange(first(r), last(r))
    end
    return res
end

@inline function indexcompare(rs::SubUnitRangesSortedSet{Ti,P,Tidx}, i, j) where {Ti,P<:UnitRangesSortedVector,Tidx}
    if i < j
        return -1
    elseif i > j
        return 1
    else
        return 0
    end
end
@inline indexcompare(rs::SubUnitRangesSortedSet{Ti,P,Tidx}, i, j) where {Ti,P<:UnitRangesSortedSet,Tidx} =
    compare(rs.parent.ranges, i, j)


@inline function findrange(rs::AbstractUnitRangesSortedSet{Ti}, i) where Ti
    if in(i, rs)
        return UnitRange{Ti}(first(i), last(i))
    else
        return nothing
    end
end

@inline getindex_tuple(rs::UnitRangesSortedVector, i) = (rs.rstarts[i], rs.rstops[i])
@inline getindex_tuple(rs::UnitRangesSortedSet, i::DataStructures.Tokens.IntSemiToken) = tuple(deref((rs.ranges, i))...)
@inline getindex_tuple(rs::Sub0UnitRangesSortedSet, i) = (rs.start, rs.stop)
@inline getindex_tuple(rs::Sub1UnitRangesSortedSet, i) = (r = getindex(rs, i); tuple(first(r), last(r)))
@inline getindex_tuple(rs::SubUnitRangesSortedSet, i) = (r = getindex(rs, i); tuple(first(r), last(r)))
@inline Base.getindex(ur::URSSUnitRange, i) = getindex(ur.parent, i)
@inline Base.getindex(rs::UnitRangesSortedVector, i) = UnitRange(rs.rstarts[i], rs.rstops[i])
@inline Base.getindex(rs::UnitRangesSortedSet, i::DataStructures.Tokens.IntSemiToken) = UnitRange(deref((rs.ranges, i))...)
@inline Base.getindex(rs::Sub0UnitRangesSortedSet, i) = rs.data
@inline function Base.getindex(rs::Sub1UnitRangesSortedSet, i)
    if i == rs.firstindex
        return rs.parent
    else
        return UnitRange(rs.stop, rs.start)
    end
end
@inline function Base.getindex(rs::SubUnitRangesSortedSet, i)
    if indexcompare(rs, i, rs.firstindexnext) == -1
        return getindexhelper_firstrange(rs, i)
    elseif indexcompare(rs, i, rs.lastindex) == -1
        return getindex(rs.parent, i)
    else
        return getindexhelper_lastrange(rs, i)
    end
end

function getindexhelper_firstrange(rs::SubUnitRangesSortedSet, i)
    if indexcompare(rs, i, rs.firstindex) == -1
        return UnitRange(rs.start + 1, rs.start) # zero range
    else
        r = getindex(rs.parent, i)
        if first(r) < rs.start
            r = UnitRange(rs.start, last(r))
        end
        if rs.stop < last(r)
            r = UnitRange(first(r), rs.stop)
        end
        return r
    end
end
function getindexhelper_lastrange(rs::SubUnitRangesSortedSet, i)
    if indexcompare(rs, rs.lastindex, i) == -1
        return UnitRange(rs.stop, rs.stop - 1) # zero range
    else
        r = getindex(rs.parent, i)
        if rs.stop < last(r)
            r = UnitRange(first(r), rs.stop)
        end
        if first(r) < rs.start
            r = UnitRange(rs.start, last(r))
        end
        return r
    end
end


@inline getindex_rangestart(rs::UnitRangesSortedVector, i) = rs.rstarts[i]
@inline getindex_rangestop(rs::UnitRangesSortedVector, i) = rs.rstops[i]
@inline getindex_rangestart(rs::UnitRangesSortedSet, i) = deref_key((rs.ranges, i))
@inline getindex_rangestop(rs::UnitRangesSortedSet, i) = deref_value((rs.ranges, i))

@inline Base.firstindex(rs::UnitRangesSortedVector) = firstindex(rs.rstarts)
@inline Base.firstindex(rs::UnitRangesSortedSet) = startof(rs.ranges)
@inline Base.firstindex(rs::Union{Sub0UnitRangesSortedSet,Sub1UnitRangesSortedSet,SubUnitRangesSortedSet}) = rs.firstindex
@inline Base.lastindex(rs::UnitRangesSortedVector) = lastindex(rs.rstarts)
@inline Base.lastindex(rs::AbstractUnitRangesSortedSet) = lastindex(rs.ranges)
@inline Base.lastindex(rs::Union{Sub0UnitRangesSortedSet,Sub1UnitRangesSortedSet,SubUnitRangesSortedSet}) = rs.lastindex
@inline Base.first(rs::UnitRangesSortedVector) = UnitRange(rs.rstarts[1], rs.rstops[1])
@inline Base.first(rs::UnitRangesSortedSet) = UnitRange(deref((rs.ranges, startof(rs.ranges)))...)
@inline Base.first(rs::Sub0UnitRangesSortedSet) = rs.parent
@inline Base.first(rs::Sub1UnitRangesSortedSet) = rs.parent
@inline function Base.first(rs::SubUnitRangesSortedSet{Ti,P,Tidx}) where {Ti,P,Tidx}
    if length(rs) != 0
        return getindex(rs, firstindex(rs))
    elseif rs.start != typemin(Ti)
        return UnitRange(rs.start, rs.start - 1)
    else
        return UnitRange(rs.start + 1, rs.start)
    end
end
@inline Base.last(rs::UnitRangesSortedVector) = UnitRange(rs.rstarts[end], rs.rstops[end])
@inline Base.last(rs::UnitRangesSortedSet) = UnitRange(deref((rs.ranges, lastindex(rs.ranges)))...)
@inline Base.last(rs::Sub0UnitRangesSortedSet) = rs.parent
@inline Base.last(rs::Sub1UnitRangesSortedSet) = rs.parent
@inline function Base.last(rs::SubUnitRangesSortedSet{Ti,P,Tidx}) where {Ti,P,Tidx}
    if length(rs) != 0
        getindex(rs, lastindex(rs))
    elseif rs.stop != typemax(Ti)
        return UnitRange(rs.stop + 1, rs.stop)
    else
        return UnitRange(rs.stop, rs.stop - 1)
    end
end

@inline beforestartindex(rs::UnitRangesSortedVector) = firstindex(rs.rstarts) - 1
@inline beforestartindex(rs::UnitRangesSortedSet) = beforestartsemitoken(rs.ranges)
@inline beforestartindex(rs::Union{Sub0UnitRangesSortedSet,Sub1UnitRangesSortedSet,SubUnitRangesSortedSet}) = rs.beforestartindex
@inline pastendindex(rs::UnitRangesSortedVector) = lastindex(rs.rstarts) + 1
@inline pastendindex(rs::UnitRangesSortedSet) = pastendsemitoken(rs.ranges)
@inline pastendindex(rs::Union{Sub0UnitRangesSortedSet,Sub1UnitRangesSortedSet,SubUnitRangesSortedSet}) = rs.pastendindex

@inline DataStructures.advance(rs::UnitRangesSortedVector, state) = state + 1
@inline DataStructures.advance(rs::UnitRangesSortedSet, state) = advance((rs.ranges, state))
@inline DataStructures.advance(rs::Union{Sub0UnitRangesSortedSet,Sub1UnitRangesSortedSet,SubUnitRangesSortedSet}, state) =
    advance(rs.parent, state)
@inline DataStructures.regress(rs::UnitRangesSortedVector, state) = state - 1
@inline DataStructures.regress(rs::UnitRangesSortedSet, state) = regress((rs.ranges, state))
@inline DataStructures.regress(rs::Union{Sub0UnitRangesSortedSet,Sub1UnitRangesSortedSet,SubUnitRangesSortedSet}, state) =
    regress(rs.parent, state)

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


#@inline searchsortedlastrange(rs::UnitRangesSortedVector, k) = searchsortedlast(rs.rstarts, k; lt=<)
"Returns index of range in which, or after, `k` is placed."
@inline searchsortedlastrange(rs::UnitRangesSortedVector, k) = bisectionsearchlast(rs.rstarts, k)
@inline searchsortedlastrange(rs::UnitRangesSortedSet, k) = searchsortedlast(rs.ranges, k)
@inline searchsortedlastrange(rs::Sub0UnitRangesSortedSet, k) = rs.beforestartindex
@inline function searchsortedlastrange(rs::Sub1UnitRangesSortedSet, k)
    if k >= rs.start
        return rs.firstindex
    else
        return rs.beforestartindex
    end
end
@inline searchsortedlastrange(rs::SubUnitRangesSortedSet{Ti,P,Tidx}, k) where {Ti,P<:UnitRangesSortedVector,Tidx} =
    bisectionsearchlast(rs.parent.rstarts, k, rs.firstindex, rs.lastindex)
@inline function searchsortedlastrange(rs::SubUnitRangesSortedSet{Ti,P,Tidx}, k) where {Ti,P<:UnitRangesSortedSet,Tidx}
    st = searchsortedlast(rs.parent.ranges, k)
    if compare(rs.parent.ranges, st, rs.beforestartindex) == 1 && compare(rs.parent.ranges, st, rs.pastendindex) == -1
        return st
    else
        return beforestartindex(rs)
    end
end

"Returns index of range in which, or before, `k` is placed."
@inline searchsortedfirstrange(rs::UnitRangesSortedVector, k) = bisectionsearchfirst(rs.rstops, k)
@inline function searchsortedfirstrange(rs::UnitRangesSortedSet, k)
    st = searchsortedlastrange(rs, k)
    if st != beforestartindex(rs) && in(k, getindex(rs, st))
        return st
    else
        return advance(rs, st)
    end
end
@inline searchsortedfirstrange(rs::Sub0UnitRangesSortedSet, k) = rs.pastendindex
@inline function searchsortedfirstrange(rs::Sub1UnitRangesSortedSet, k)
    if k <= rs.stop
        return rs.lastindex
    else
        return rs.pastendindex
    end
end
@inline searchsortedfirstrange(rs::SubUnitRangesSortedSet{Ti,P,Tidx}, k) where {Ti,P<:UnitRangesSortedVector,Tidx} =
    bisectionsearchfirst(rs.parent.rstops, k, rs.firstindex, rs.lastindex)
@inline function searchsortedfirstrange(rs::SubUnitRangesSortedSet{Ti,P,Tidx}, k) where {Ti,P<:UnitRangesSortedSet,Tidx}
    st = searchsortedlastrange(rs, k)
    if st != beforestartindex(rs) && in(k, getindex(rs, st))
        return st
    else
        return advance(rs, st)
    end
end

"Returns range of `rs` indexes which coincide or concluded in `I` range."
@inline searchsortedrange(rs::UnitRangesSortedVector, I::AbstractRange) =
    UnitRange(searchsortedfirstrange(rs, first(I)), searchsortedlastrange(rs, last(I)))
@inline searchsortedrange(rs::UnitRangesSortedSet, I::AbstractRange) =
    URSSUnitRange(rs, searchsortedfirstrange(rs, first(I)), searchsortedlastrange(rs, last(I)))
@inline searchsortedrange(rs::Sub0UnitRangesSortedSet{Ti,P,Tidx}, I::AbstractRange) where {Ti,P<:UnitRangesSortedVector,Tidx} =
    UnitRange(searchsortedfirstrange(rs, first(I)), searchsortedlastrange(rs, last(I)))
@inline searchsortedrange(rs::Sub0UnitRangesSortedSet{Ti,P,Tidx}, I::AbstractRange) where {Ti,P<:UnitRangesSortedSet,Tidx} =
    URSSUnitRange(rs, searchsortedfirstrange(rs, first(I)), searchsortedlastrange(rs, last(I)))
@inline searchsortedrange(rs::Sub1UnitRangesSortedSet{Ti,P,Tidx}, I::AbstractRange) where {Ti,P<:UnitRangesSortedVector,Tidx} =
    UnitRange(searchsortedfirstrange(rs, first(I)), searchsortedlastrange(rs, last(I)))
@inline searchsortedrange(rs::Sub1UnitRangesSortedSet{Ti,P,Tidx}, I::AbstractRange) where {Ti,P<:UnitRangesSortedSet,Tidx} =
    URSSUnitRange(rs, searchsortedfirstrange(rs, first(I)), searchsortedlastrange(rs, last(I)))
@inline searchsortedrange(rs::SubUnitRangesSortedSet{Ti,P,Tidx}, I::AbstractRange) where {Ti,P<:UnitRangesSortedVector,Tidx} =
    UnitRange(searchsortedfirstrange(rs, first(I)), searchsortedlastrange(rs, last(I)))
@inline searchsortedrange(rs::SubUnitRangesSortedSet{Ti,P,Tidx}, I::AbstractRange) where {Ti,P<:UnitRangesSortedSet,Tidx} =
    URSSUnitRange(rs, searchsortedfirstrange(rs, first(I)), searchsortedlastrange(rs, last(I)))

"Returns indexes of range in `rs` in which `k` may be inserted. Or negative range in the case of `k` is
 between `rs` ranges, and indices of resulted range is the indexes of that neighbors."
@inline function searchsortedrange(rs::UnitRangesSortedVector{Ti}, k) where Ti
    st = searchsortedlastrange(rs, k)
    if st != beforestartindex(rs) && in(k, getindex(rs, st))
        return UnitRange(st, st)
    else
        return UnitRange(advance(rs, st), st)
    end
end
@inline function searchsortedrange(rs::UnitRangesSortedSet{Ti}, k) where Ti
    st = searchsortedlastrange(rs, k)
    if st != beforestartindex(rs) && in(k, getindex(rs, st))
        return URSSUnitRange(rs, st, st)
    else
        return URSSUnitRange(rs, advance(rs, st), st)
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

#
# Assignments
#

@inline function Base.in(II::AbstractRange, rs::AbstractUnitRangesSortedSet{Ti}) where Ti
    I = UnitRange{Ti}(first(II), last(II))
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
@inline Base.in(II::AbstractRange, su::Sub1UnitRangesSortedSet{Ti,P,Tidx}) where {Ti,P,Tidx} =
    issubset(UnitRange{Ti}(first(II), last(II)), su.data)
@inline Base.in(key, su::Sub1UnitRangesSortedSet{Ti,P,Tidx}) where {Ti,P,Tidx} = in(Ti(key), su.data)

@inline function Base.in(II::AbstractRange, su::SubUnitRangesSortedSet{Ti,P,Tidx}) where {Ti,P,Tidx}
    I = UnitRange{Ti}(first(II), last(II))
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

@inline function Base.in(key, su::SubUnitRangesSortedSet{Ti,P,Tidx}) where {Ti,P,Tidx}
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

function Base.push!(rs::AbstractUnitRangesSortedSet, II::AbstractUnitRangesSortedSet)
    for r in II
        push!(rs, r)
    end
    rs
end

function Base.push!(rs::UnitRangesSortedVector{Ti}, II::AbstractRange) where Ti
    I = UnitRange{Ti}(first(II), last(II))

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
    @boundscheck @assert Iposition == searchsortedrange(rs, last(I)) "FIXME: I Am an error in program logic."

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
        throw(AssertionError("FIXME: I Am an error in program logic."))
    end

    return rs
end

function Base.push!(rs::UnitRangesSortedSet{Ti}, II::AbstractRange) where Ti
    I = UnitRange{Ti}(first(II), last(II))

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
    @boundscheck @assert Iposition == searchsortedrange(rs, last(I)) "FIXME: I Am an error in program logic."

    Iposition_left = last(Iposition)
    Iposition_right = first(Iposition)
    Iposition_left_range_stop = Iposition_left != beforestartindex(rs) ? getindex_rangestop(rs, Iposition_left) :
                                                                         typemin(Ti)
    Iposition_right_range_start = Iposition_right != pastendindex(rs) ? getindex_rangestart(rs, Iposition_right) :
                                                                        typemax(Ti)

    # `I` is adjoined in both sides with `rs` ranges, thus join them all
    if Iposition_left != beforestartindex(rs) && Iposition_left_range_stop + 1 == first(I) &&
       Iposition_right != pastendindex(rs) && Iposition_right_range_start - 1 == last(I)
        rs.ranges[Iposition_left] = getindex_rangestop(rs, Iposition_right)
        delete!((rs.ranges, Iposition_right))

    # `I` is adjoin with `rs` range on left side, thus append to left range
    elseif Iposition_left != beforestartindex(rs) && Iposition_left_range_stop + 1 == first(I) &&
           Iposition_right != pastendindex(rs) && Iposition_right_range_start - 1 != last(I)
        rs.ranges[Iposition_left] = last(I)

    # `I` is adjoin with `rs` range on right side, thus prepend to right range
    elseif Iposition_left != beforestartindex(rs) && Iposition_left_range_stop + 1 != first(I) &&
           Iposition_right != pastendindex(rs) && Iposition_right_range_start - 1 == last(I)
        rs.ranges[first(I)] = getindex_rangestop(rs, Iposition_right)
        delete!((rs.ranges, Iposition_right))

    # `I` is separate from both sides, insert it
    elseif Iposition_left != beforestartindex(rs) && Iposition_left_range_stop + 1 != first(I) &&
           Iposition_right != pastendindex(rs) && Iposition_right_range_start - 1 != last(I)
        rs.ranges[first(I)] = last(I)

    # `I` is first range in `rs` and adjoin with range on right side, thus prepend `I` to right range
    elseif Iposition_left == beforestartindex(rs) &&
           Iposition_right != pastendindex(rs) && Iposition_right_range_start - 1 == last(I)
        rs.ranges[first(I)] = getindex_rangestop(rs, Iposition_right)
        delete!((rs.ranges, Iposition_right))

    # `I` is separate first range in `rs`, insert it
    elseif Iposition_left == beforestartindex(rs) &&
           Iposition_right != pastendindex(rs) && Iposition_right_range_start - 1 != last(I)
        rs.ranges[first(I)] = last(I)

    # `I` is last range in `rs` and adjoin with range on left side, thus append `I` to left range
    elseif Iposition_left != beforestartindex(rs) && Iposition_left_range_stop + 1 == first(I) &&
           Iposition_right == pastendindex(rs)
        rs.ranges[Iposition_left] = last(I)

    # `I` is separate last range in `rs`, insert it
    elseif Iposition_left != beforestartindex(rs) && Iposition_left_range_stop + 1 != first(I) &&
           Iposition_right == pastendindex(rs)
        rs.ranges[first(I)] = last(I)

    # `rs` is empty
    elseif Iposition_left == beforestartindex(rs) && Iposition_right == pastendindex(rs)
        rs.ranges[first(I)] = last(I)

    else
        throw(AssertionError("FIXME: I Am an error in program logic."))
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
        throw(KeyError(k))
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


function Base.push!(rs::UnitRangesSortedVector{Ti}, key) where {Ti}
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


function Base.push!(rs::UnitRangesSortedSet{Ti}, key) where {Ti}
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


function Base.delete!(rs::AbstractUnitRangesSortedSet, II::AbstractUnitRangesSortedSet)
    for r in II
        delete!(rs, r)
    end
    rs
end

function Base.delete!(rs::UnitRangesSortedVector{Ti}, II::AbstractRange) where {Ti}
    length(II) == 0 && return rs

    I = UnitRange{Ti}(first(II), last(II))

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

    Iposition_range_start = getindex_rangestart(rs, first(Iposition))
    Iposition_range_stop = getindex_rangestop(rs, last(Iposition))

    # `I` intersects only one range in `rs`
    if length(Iposition) == 1

        # inside one range, thus split it
        if Iposition_range_start < first(I) && last(I) < Iposition_range_stop
            insert!(rs.rstarts, advance(rs, first(Iposition)), last(I) + 1)
            insert!(rs.rstops, first(Iposition), first(I) - 1)

        # `I` intersects with or inside in one range from left side, thus shrink from left side
        elseif Iposition_range_start >= first(I) && last(I) < Iposition_range_stop
            rs.rstarts[first(Iposition)] = last(I) + 1

        # `I` intersects with or inside in one range from right side, thus shrink from right side
        # inside one range, ended to left side, thus shrink from right side
        elseif Iposition_range_start < first(I) && last(I) >= Iposition_range_stop
            rs.rstops[first(Iposition)] = first(I) - 1

        else
            throw(AssertionError("FIXME: I Am an error in program logic."))
        end

    # All remaining cases are with two ranges from `rs`,
    # and both side ranges only partly intersects with `I`.
    elseif length(Iposition) == 2

        # shrink second (right) range from the left side
        rs.rstarts[last(Iposition)] = last(I) + 1

        # shrink first (left) range from the right side
        rs.rstops[first(Iposition)] = first(I) - 1

    else
        throw(AssertionError("FIXME: I Am an error in program logic."))
    end

    return rs
end


function Base.delete!(rs::UnitRangesSortedSet{Ti}, II::AbstractRange) where {Ti}
    length(II) == 0 && return rs

    I = UnitRange{Ti}(first(II), last(II))

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

    Iposition_range_start = getindex_rangestart(rs, first(Iposition))
    Iposition_range_stop = getindex_rangestop(rs, last(Iposition))

    # `I` intersects only one range in `rs`
    if length(Iposition) == 1

        # inside one range, thus split it
        if Iposition_range_start < first(I) && last(I) < Iposition_range_stop
            rs.ranges[last(I) + 1] = Iposition_range_stop
            rs.ranges[first(Iposition)] = first(I) - 1

        # `I` intersects with or inside in one range from left side, thus shrink from left side
        elseif Iposition_range_start >= first(I) && last(I) < Iposition_range_stop
            rs.ranges[last(I) + 1] = Iposition_range_stop
            delete!((rs.ranges, first(Iposition)))

        # `I` intersects with or inside in one range from right side, thus shrink from right side
        # inside one range, ended to left side, thus shrink from right side
        elseif Iposition_range_start < first(I) && last(I) >= Iposition_range_stop
            rs.ranges[first(Iposition)] = first(I) - 1

        else
            throw(AssertionError("FIXME: I Am an error in program logic."))
        end

    # All remaining cases are with two ranges from `rs`,
    # and both side ranges only partly intersects with `I`.
    elseif length(Iposition) == 2

        # shrink second (right) range from the left side
        rs.ranges[last(I) + 1] = Iposition_range_stop
        delete!((rs.ranges, last(Iposition)))


        # shrink first (left) range from the right side
        rs.ranges[first(Iposition)] = first(I) - 1

    else
        throw(AssertionError("FIXME: I Am an error in program logic."))
    end

    return rs
end


function Base.delete!(rs::UnitRangesSortedVector{Ti}, key) where {Ti}
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


function Base.delete!(rs::UnitRangesSortedSet{Ti}, key) where {Ti}
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

function Base.pop!(rs::AbstractUnitRangesSortedSet, k)
    if (r = findrange(rs, k)) !== nothing
        delete!(rs, k)
        return r
    else
        throw(KeyError(k))
    end
end
function Base.pop!(rs::AbstractUnitRangesSortedSet, k, default)
    if (r = findrange(rs, k)) !== nothing
        delete!(rs, k)
        return r
    else
        return default
    end
end


# TODO: setdiff, symdiff, eachindex (SomeEachindexIterator), map
# setdiff already exist in julia/base/abstractset.jl


Base.empty(rs::T) where {T<:AbstractUnitRangesSortedSet} = T()

function Base.empty!(rs::UnitRangesSortedVector)
    empty!(rs.rstarts)
    empty!(rs.rstops)
    rs.lastusedrangeindex = beforestartindex(rs)
    rs
end
function Base.empty!(rs::AbstractUnitRangesSortedSet)
    empty!(rs.ranges)
    rs.lastusedrangeindex = beforestartindex(rs)
    rs
end

Base.copy(rs::T) where {T<:UnitRangesSortedVector} = T(rs.lastusedrangeindex, copy(rs.rstarts), copy(rs.rstops))
Base.copy(rs::T) where {T<:UnitRangesSortedSet} = T(rs.lastusedrangeindex, packcopy(rs.ranges))

#Base.deepcopy(rs::AbstractUnitRangesSortedSet) = copy(rs)



Base.union(rs::AbstractUnitRangesSortedSet, rss...) = union!(copy(rs), rss...)
Base.union(rs::AbstractUnitRangesSortedSet, rs2) = union!(copy(rs), rs2)

@inline Base.union!(rs::AbstractUnitRangesSortedSet, rss...) = union!(union!(rs, rss[1]), Base.tail(rss)...)
function Base.union!(rs::AbstractUnitRangesSortedSet, rs2)
    issubset(rs2, rs) && return rs
    for r in rs2
        push!(rs, r)
    end
    return rs
end

@inline Base.union!(rs::AbstractUnitRangesSortedSet, r::AbstractRange) = push!(rs, r)

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
function Base.issubset(rs1::AbstractUnitRangesSortedSet, rs2::Union{AbstractSet{T},AbstractVector{T},Tuple{T}}) where {T<:AbstractRange}
    for r1 in rs1
        findfirst(s->issubset(r1, s), rs2) != nothing || return false
    end
    return true
end


function Base.filter(pred::Function, rs::AbstractUnitRangesSortedSet)
    res = empty(rs)
    for r in rs
        pred(r) && push!(res, r)
    end
    res
end
function Base.filter!(pred::Function, rs::AbstractUnitRangesSortedSet)
    for r in Iterators.reverse(rs)
        pred(r) || delete!(rs, r)
    end
    rs
end

function Base.intersect!(rs::AbstractUnitRangesSortedSet{Ti}, II::AbstractRange) where Ti
    length(II) == 0 && return rs
    I = UnitRange{Ti}(first(II), last(II))

    delete!(rs, getindex_rangestart(rs, firstindex(rs)):first(I) - 1)
    delete!(rs, last(I) + 1:getindex_rangestop(rs, lastindex(rs)))

    rs
end

function Base.intersect!(rs1::AbstractUnitRangesSortedSet, rs2::AbstractUnitRangesSortedSet)
    length(rs2) == 0 && return rs1

    intersect!(rs1, getindex_rangestart(rs2, firstindex(rs2)):getindex_rangestop(rs2, lastindex(rs2)))

    length(rs2) == 1 && return rs1

    rprev, rs2tail = Iterators.peel(rs2)
    for r in rs2tail
        delete!(rs1, last(rprev)+1:first(r)-1)
        rprev = r
    end

    rs1
end

Base.intersect(rs1::AbstractUnitRangesSortedSet, rs2::AbstractUnitRangesSortedSet) = intersect!(copy(rs1), rs2)


#
#  Aux functions
#

Base.show(io::IO, ur::URSSUnitRange) = print(io, repr(first(ur)), ':', length(ur) == 0 ? "0" : "1" , ':', repr(last(ur)))

function Base.show(io::IO, ::MIME"text/plain", x::URSSUnitRange)
    len = length(x)
    print(io, len, "-element range ", first(x), ":", last(x), " in ", typeof(x.parent))
    if len != 0
        println(io, " containing:")
        show(IOContext(io, :typeinfo => eltype(x)), x)
    end
end

function Base.show(io::IOContext, x::URSSUnitRange)
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


function Base.show(io::IO, x::AbstractUnitRangesSortedSet)
    print(io, typeof(x), "(")
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
function Base.show(io::IO, ::MIME"text/plain", x::AbstractUnitRangesSortedSet)
    len = length(x)
    print(io, typeof(x), "()")
    if len != 0
        println(io, ":")
        show(IOContext(io, :typeinfo => eltype(x)), x)
    end
end

function Base.show(io::IOContext, x::AbstractUnitRangesSortedSet)
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
