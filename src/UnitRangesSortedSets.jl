
module UnitRangesSortedSets
export AbstractUnitRangesSortedSet, UnitRangesSortedVector, UnitRangesSortedSet, SubUnitRangesSortedSet, URSSUnitRange
export iteratenzpairs, iteratenzpairsview, iteratenzvalues, iteratenzvaluesview, iteratenzindices
export testfun_create, testfun_createSV, testfun_createVL, testfun_create_seq, testfun_create_dense, testfun_delete!, testfun_in, testfun_in_outer, testfun_in_rand, testfun_in_seq, testfun_nzgetindex, testfun_setindex!
export searchsortedrange, searchsortedfirstrange, searchsortedlastrange


import Base: ForwardOrdering, Forward
const FOrd = ForwardOrdering

using DocStringExtensions
using DataStructures
using IterTools
using Setfield
#using SparseArrays
#using StaticArrays
#import SparseArrays: indtype, nonzeroinds, nonzeros
using Random


abstract type AbstractUnitRangesSortedSet{Ti} <: AbstractSet{Ti} end

struct SubUnitRangesSortedSet{Ti,P} <: AbstractUnitRangesSortedSet{Ti}
    parent::P
    start::Ti
    stop::Ti
    startindex::Int
    stopindex::Int
end

function Base.view(rs::Tp, I::AbstractRange) where {Tp<:AbstractUnitRangesSortedSet{Ti}} where Ti
    r = searchsortedrange(rs, I)
    return SubUnitRangesSortedSet{Ti,Tp}(rs, Ti(first(I)), Ti(last(I)), first(r), last(r))
end

"""
Inserting zero, or negative length ranges does nothing.
$(TYPEDEF)
Mutable struct fields:
$(TYPEDFIELDS)
"""
mutable struct UnitRangesSortedVector{Ti} <: AbstractUnitRangesSortedSet{Ti}
    "Index of last used chunk"
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
    "Index of last used chunk"
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

@inline Base.:(==)(l::URSSUnitRange{Trs,Tidx}, r::URSSUnitRange{Trs,Tidx}) where {Trs<:UnitRangesSortedVector,Tidx} =
    l.parent === r.parent &&
    l.start == r.start &&
    l.stop == r.stop
@inline Base.:(==)(l::URSSUnitRange{Trs,Tidx}, r::URSSUnitRange{Trs,Tidx}) where {Trs<:UnitRangesSortedSet,Tidx} =
    l.parent === r.parent &&
    compare(l.parent.ranges, l.start, r.start) == 0 &&
    compare(l.parent.ranges, l.stop, r.stop) == 0

@inline Base.firstindex(ur::URSSUnitRange) = ur.start
@inline Base.lastindex(ur::URSSUnitRange) = ur.stop
@inline DataStructures.advance(ur::URSSUnitRange, i::DataStructures.Tokens.IntSemiToken) = advance((ur.parent, i))
@inline DataStructures.regress(ur::URSSUnitRange, i::DataStructures.Tokens.IntSemiToken) = regress((ur.parent, i))
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
#    map(s->get_range(ur.parent, s), ur.start:ur.stop)
#@inline Base.collect(ur::URSSUnitRange{Trs,Tidx}) where {Trs<:UnitRangesSortedSet,Tidx} =
#    map(s->get_range(ur.parent, s), onlysemitokens(inclusive(ur.parent.ranges, ur.start, ur.stop)))

Base.length(rs::UnitRangesSortedVector) = length(rs.rstarts)
Base.length(rs::UnitRangesSortedSet) = length(rs.ranges)
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

@inline get_range_length(rs::AbstractUnitRangesSortedSet, i) = length(get_range(rs, i))
@inline get_range_indices(rs::UnitRangesSortedVector, i) = (rs.rstarts[i], rs.rstops[i])
@inline get_range_indices(rs::UnitRangesSortedSet, i::DataStructures.Tokens.IntSemiToken) = tuple(deref((rs.ranges, i))...)
@inline get_range(rs::UnitRangesSortedVector, i) = UnitRange(rs.rstarts[i], rs.rstops[i])
@inline get_range(rs::UnitRangesSortedSet, i::DataStructures.Tokens.IntSemiToken) = UnitRange(deref((rs.ranges, i))...)
@inline get_range(ur::URSSUnitRange, i) = get_range(ur.parent, i)
function get_range(rs::SubArray{<:Any,<:Any,<:T}, i) where {T<:AbstractUnitRangesSortedSet}
    idx1 = first(rs.indices[1])
    key = rs.parent.ranges[i]
    len = length(rs.parent.ranges[i])
    if key <= idx1 < key + len
        return @view(rs.parent.ranges[i][idx1:end])
    elseif key <= last(rs.indices[1]) < key + len
        return view(rs.parent.ranges[i], 1:(last(rs.indices[1])-key+1))
    else
        return @view(rs.parent.ranges[i][1:end])
    end
end
function get_range(rs::SubArray{<:Any,<:Any,<:T}, i) where {T<:AbstractUnitRangesSortedSet}
    idx1 = first(rs.indices[1])
    key, chunk = deref((rs.parent.ranges, i))
    len = length(chunk)
    if key <= idx1 < key + len
        return @view(chunk[idx1:end])
    elseif key <= last(rs.indices[1]) < key + len
        return view(chunk, 1:(last(rs.indices[1])-key+1))
    else
        return @view(chunk[1:end])
    end
end
@inline get_range_start(rs::UnitRangesSortedVector, i) = rs.rstarts[i]
@inline get_range_stop(rs::UnitRangesSortedVector, i) = rs.rstops[i]
@inline get_range_start(rs::UnitRangesSortedSet, i) = deref_key((rs.ranges, i))
@inline get_range_stop(rs::UnitRangesSortedSet, i) = deref_value((rs.ranges, i))
@inline function get_range_start(rs::SubArray{<:Any,<:Any,<:T}, i) where {T<:AbstractUnitRangesSortedSet}
    if rs.parent.ranges[i] <= first(rs.indices[1]) < rs.parent.ranges[i] + length(rs.parent.ranges[i])
        return first(rs.indices[1])
    else
        return rs.parent.ranges[i]
    end
end
@inline function get_range_start(rs::SubArray{<:Any,<:Any,<:T}, i) where {T<:AbstractUnitRangesSortedSet}
    key, chunk = deref((rs.parent.ranges, i))
    len = length(chunk)
    if key <= first(rs.indices[1]) < key + len
        return first(rs.indices[1])
    else
        return key
    end
end

@inline Base.firstindex(rs::UnitRangesSortedVector) = firstindex(rs.rstarts)
@inline Base.firstindex(rs::UnitRangesSortedSet) = startof(rs.ranges)
@inline Base.lastindex(rs::UnitRangesSortedVector) = lastindex(rs.rstarts)
@inline Base.lastindex(rs::AbstractUnitRangesSortedSet) = lastindex(rs.ranges)
@inline Base.first(rs::UnitRangesSortedVector) = UnitRange(rs.rstarts[1], rs.rstops[1])
@inline Base.first(rs::UnitRangesSortedSet) = UnitRange(deref((rs.ranges, startof(rs.ranges)))...)
@inline Base.last(rs::UnitRangesSortedVector) = UnitRange(rs.rstarts[end], rs.rstops[end])
@inline Base.last(rs::UnitRangesSortedSet) = UnitRange(deref((rs.ranges, lastindex(rs.ranges)))...)

@inline beforestartindex(rs::UnitRangesSortedVector) = firstindex(rs.rstarts) - 1
@inline beforestartindex(rs::UnitRangesSortedSet) = beforestartsemitoken(rs.ranges)
@inline pastendindex(rs::UnitRangesSortedVector) = lastindex(rs.rstarts) + 1
@inline pastendindex(rs::UnitRangesSortedSet) = pastendsemitoken(rs.ranges)

@inline DataStructures.advance(rs::UnitRangesSortedVector, state) = state + 1
@inline DataStructures.advance(rs::UnitRangesSortedSet, state) = advance((rs.ranges, state))
@inline DataStructures.regress(rs::UnitRangesSortedVector, state) = state - 1
@inline DataStructures.regress(rs::UnitRangesSortedSet, state) = regress((rs.ranges, state))

# Derived from julia/base/sort.jl, for increasing sorted vectors.
function bisectionsearchlast(V::AbstractVector, i)
    lo = 0
    hi = length(V) + 1
    @inbounds while lo < hi - 1
        m = lo + ((hi - lo) >>> 0x01)
        if i < V[m]
            hi = m
        else
            lo = m
        end
    end
    return lo
end
function bisectionsearchlast(V::AbstractVector, i, lo::T, hi::T) where T
    u = T(1)
    lo = lo - u
    hi = hi + u
    @inbounds while lo < hi - u
        m = lo + ((hi - lo) >>> 0x01)
        if i < V[m]
            hi = m
        else
            lo = m
        end
    end
    return lo
end

function bisectionsearchfirst(V::AbstractVector, i)
    lo = 0
    hi = length(V) + 1
    @inbounds while lo < hi - 1
        m = lo + ((hi - lo) >>> 0x01)
        if V[m] < i
            lo = m
        else
            hi = m
        end
    end
    return hi
end
function bisectionsearchfirst(V::AbstractVector, i, lo::T, hi::T) where T
    u = T(1)
    lo = lo - u
    hi = hi + u
    @inbounds while lo < hi - u
        m = lo + ((hi - lo) >>> 0x01)
        if V[m] < i
            lo = m
        else
            hi = m
        end
    end
    return hi
end


#@inline searchsortedlastrange(rs::UnitRangesSortedVector, i) = searchsortedlast(rs.rstarts, i; lt=<)
"Returns index of range in which, or after, `i` is placed."
@inline searchsortedlastrange(rs::UnitRangesSortedVector, i) = bisectionsearchlast(rs.rstarts, i)
@inline searchsortedlastrange(rs::UnitRangesSortedSet, i) = searchsortedlast(rs.ranges, i)
"Returns index of range in which, or before, `i` is placed."
@inline searchsortedfirstrange(rs::UnitRangesSortedVector, i) = bisectionsearchfirst(rs.rstops, i)
@inline function searchsortedfirstrange(rs::UnitRangesSortedSet, i)
    st = searchsortedlastrange(rs, i)
    if st != beforestartindex(rs) && in(i, get_range(rs, st))
        return st
    else
        return advance(rs, st)
    end
end

"Returns range of `rs` indexes which coincide or concluded in `I` range."
@inline searchsortedrange(rs::UnitRangesSortedVector{Ti}, I::AbstractRange) where Ti =
    UnitRange(searchsortedfirstrange(rs, first(I)), searchsortedlastrange(rs, last(I)))
@inline searchsortedrange(rs::UnitRangesSortedSet{Ti}, I::AbstractRange) where Ti =
    URSSUnitRange(rs, searchsortedfirstrange(rs, first(I)), searchsortedlastrange(rs, last(I)))

"Returns indexes of range in `rs` in which `i` may be inserted. Or negative range in the case of `i` is
 between `rs` ranges, and indices of resulted range is the indexes of that neighbors."
@inline function searchsortedrange(rs::UnitRangesSortedVector{Ti}, i) where Ti
    st = searchsortedlastrange(rs, i)
    if st != beforestartindex(rs) && in(i, get_range(rs, st))
        return UnitRange(st, st)
    else
        return UnitRange(advance(rs, st), st)
    end
end
@inline function searchsortedrange(rs::UnitRangesSortedSet{Ti}, i) where Ti
    st = searchsortedlastrange(rs, i)
    if st != beforestartindex(rs) && in(i, get_range(rs, st))
        return URSSUnitRange(rs, st, st)
    else
        return URSSUnitRange(rs, advance(rs, st), st)
    end
end


"Returns nzchunk which on vector index `i`, or after `i`"
@inline function searchsortedlast_nzchunk(rs::AbstractUnitRangesSortedSet, i::Integer)
    if i == 1 # most of use cases
        return nnz(rs) == 0 ? pastendindex(rs) : firstindex(rs)
    elseif nnz(rs) != 0
        st = searchsortedlastchunk(rs, i)
        if st != beforestartindex(rs)
            key = get_range_start(rs, st)
            len = get_range_length(rs, st)
            if i < key + len
                return st
            else
                return advance(rs, st)
            end
        else
            return firstindex(rs)
        end
    else
        return beforestartindex(rs)
    end
end

"Returns nzchunk which on vector index `i`, or before `i`"
@inline function searchsortedfirst_nzchunk(rs::AbstractUnitRangesSortedSet, i::Integer)
    if nnz(rs) != 0
        return searchsortedlastchunk(rs, i)
    else
        return beforestartindex(rs)
    end
end

@inline index_status(rs::UnitRangesSortedVector, st) = 0 <= st <= length(rs)+1 ? 1 : 0
@inline index_status(rs::UnitRangesSortedSet, st) = status((rs.ranges, st))

#@inline function Base.findfirst(testf::Function, rs::AbstractUnitRangesSortedSet)
#    for p in nzpairs(rs)
#        testf(last(p)) && return first(p)
#    end
#    return nothing
#end

Base.findall(testf::Function, rs::AbstractUnitRangesSortedSet) = collect(p for p in rs if testf(p))


#
#  Iterators
#


@inline function Base.iterate(rs::AbstractUnitRangesSortedSet, state = firstindex(rs))
    if state != pastendindex(rs)
        return (get_range(rs, state), advance(rs, state))
    else
        return nothing
    end
end

@inline function Base.iterate(rrs::Base.Iterators.Reverse{T}, state = lastindex(rrs.itr)) where {T<:AbstractUnitRangesSortedSet}
    if state != beforestartindex(rrs.itr)
        return (get_range(rrs.itr, state), regress(rrs.itr, state))
    else
        return nothing
    end
end

#
# Assignments
#

@inline function Base.in(II::AbstractRange, rs::AbstractUnitRangesSortedSet{Ti}) where Ti
    I = convert(UnitRange{Ti}, II)
    # fast check for cached range index
    if (st = rs.lastusedrangeindex) != beforestartindex(rs)
        if issubset(I, get_range(rs, st))
            return true
        end
    end
    # cached range index miss (or index not stored), thus try search
    Iposition = searchsortedrange(rs, I)
    rs.lastusedrangeindex = last(Iposition)
    if length(Iposition) != 1
        return false
    elseif issubset(I, get_range(rs, first(Iposition)))
        return true
    else
        return false
    end
end

@inline function Base.in(idx, rs::AbstractUnitRangesSortedSet{Ti}) where Ti
    i = Ti(idx)
    # fast check for cached range index
    if (st = rs.lastusedrangeindex) != beforestartindex(rs)
        r_start, r_stop = get_range_indices(rs, st)
        if r_start <= i <= r_stop
            return true
        end
    end
    # cached range index miss (or index not stored), thus try search
    st = searchsortedlastrange(rs, i)
    if st != beforestartindex(rs)  # the index `i` is not before the start of first range
        r_start, r_stop = get_range_indices(rs, st)
        if i <= r_stop  # is the index `i` inside of range
            rs.lastusedrangeindex = st
            return true
        end
    end
    rs.lastusedrangeindex = beforestartindex(rs)
    return false
end


#@inline Base.haskey(rs::AbstractUnitRangesSortedSet, idx) = in(idx, rs)

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
    if length(Iposition) == 1 && I == get_range(rs, first(Iposition))
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
    if Iposition_left != beforestartindex(rs) && get_range_stop(rs, Iposition_left) + 1 == first(I) &&
       Iposition_right != pastendindex(rs) && get_range_start(rs, Iposition_right) - 1 == last(I)
        rs.rstops[Iposition_left] = get_range_stop(rs, Iposition_right)
        deleteat!(rs.rstarts, Iposition_right)
        deleteat!(rs.rstops, Iposition_right)

    # `I` is adjoin with `rs` range on left side, thus append to left range
    elseif Iposition_left != beforestartindex(rs) && get_range_stop(rs, Iposition_left) + 1 == first(I) &&
           Iposition_right != pastendindex(rs) && get_range_start(rs, Iposition_right) - 1 != last(I)
        rs.rstops[Iposition_left] = last(I)

    # `I` is adjoin with `rs` range on right side, thus prepend to right range
    elseif Iposition_left != beforestartindex(rs) && get_range_stop(rs, Iposition_left) + 1 != first(I) &&
           Iposition_right != pastendindex(rs) && get_range_start(rs, Iposition_right) - 1 == last(I)
        rs.rstarts[Iposition_right] = first(I)

    # `I` is separate from both sides, insert it
    elseif Iposition_left != beforestartindex(rs) && get_range_stop(rs, Iposition_left) + 1 != first(I) &&
           Iposition_right != pastendindex(rs) && get_range_start(rs, Iposition_right) - 1 != last(I)
        insert!(rs.rstarts, Iposition_right, first(I))
        insert!(rs.rstops, Iposition_right, last(I))

    # `I` is first range in `rs` and adjoin with range on right side, thus prepend `I` to right range
    elseif Iposition_left == beforestartindex(rs) &&
           Iposition_right != pastendindex(rs) && get_range_start(rs, Iposition_right) - 1 == last(I)
        rs.rstarts[Iposition_right] = first(I)

    # `I` is separate first range in `rs`, insert it
    elseif Iposition_left == beforestartindex(rs) &&
           Iposition_right != pastendindex(rs) && get_range_start(rs, Iposition_right) - 1 != last(I)
        insert!(rs.rstarts, Iposition_right, first(I))
        insert!(rs.rstops, Iposition_right, last(I))

    # `I` is last range in `rs` and adjoin with range on left side, thus append `I` to left range
    elseif Iposition_left != beforestartindex(rs) && get_range_stop(rs, Iposition_left) + 1 == first(I) &&
           Iposition_right == pastendindex(rs)
        rs.rstops[Iposition_left] = last(I)

    # `I` is separate last range in `rs`, insert it
    elseif Iposition_left != beforestartindex(rs) && get_range_stop(rs, Iposition_left) + 1 != first(I) &&
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
    if length(Iposition) == 1 && I == get_range(rs, first(Iposition))
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
    Iposition_left_range_stop = Iposition_left != beforestartindex(rs) ? get_range_stop(rs, Iposition_left) :
                                                                         typemin(Ti)
    Iposition_right_range_start = Iposition_right != pastendindex(rs) ? get_range_start(rs, Iposition_right) :
                                                                        typemax(Ti)

    # `I` is adjoined in both sides with `rs` ranges, thus join them all
    if Iposition_left != beforestartindex(rs) && Iposition_left_range_stop + 1 == first(I) &&
       Iposition_right != pastendindex(rs) && Iposition_right_range_start - 1 == last(I)
        rs.ranges[Iposition_left] = get_range_stop(rs, Iposition_right)
        delete!((rs.ranges, Iposition_right))

    # `I` is adjoin with `rs` range on left side, thus append to left range
    elseif Iposition_left != beforestartindex(rs) && Iposition_left_range_stop + 1 == first(I) &&
           Iposition_right != pastendindex(rs) && Iposition_right_range_start - 1 != last(I)
        rs.ranges[Iposition_left] = last(I)

    # `I` is adjoin with `rs` range on right side, thus prepend to right range
    elseif Iposition_left != beforestartindex(rs) && Iposition_left_range_stop + 1 != first(I) &&
           Iposition_right != pastendindex(rs) && Iposition_right_range_start - 1 == last(I)
        rs.ranges[first(I)] = get_range_stop(rs, Iposition_right)
        delete!((rs.ranges, Iposition_right))

    # `I` is separate from both sides, insert it
    elseif Iposition_left != beforestartindex(rs) && Iposition_left_range_stop + 1 != first(I) &&
           Iposition_right != pastendindex(rs) && Iposition_right_range_start - 1 != last(I)
        rs.ranges[first(I)] = last(I)

    # `I` is first range in `rs` and adjoin with range on right side, thus prepend `I` to right range
    elseif Iposition_left == beforestartindex(rs) &&
           Iposition_right != pastendindex(rs) && Iposition_right_range_start - 1 == last(I)
        rs.ranges[first(I)] = get_range_stop(rs, Iposition_right)
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

@inline function check_exist_and_update(rs, i)

    # fast check for cached range index
    if (st = rs.lastusedrangeindex) != beforestartindex(rs)
        r_start, r_stop = get_range_indices(rs, st)
        if r_start <= i <= r_stop
            return nothing
        end
    end

    st = searchsortedlastrange(rs, i)
    #
    #sstatus = status((rs.ranges, st))
    @boundscheck if index_status(rs, st) == 0 # invalid range index
        throw(KeyError(i))
    end

    # check the index exist and update its data
    if st != beforestartindex(rs)  # the index `i` is not before the first index
        if i <= get_range_stop(rs, st)
            rs.lastusedrangeindex = st
            return nothing
        end
    end

    return st
end


function Base.push!(rs::UnitRangesSortedVector{Ti}, idx) where {Ti}
    i = Ti(idx)

    st = check_exist_and_update(rs, i)
    st == nothing && return rs

    if length(rs) == 0
        push!(rs.rstarts, i)
        push!(rs.rstops, i)
        rs.lastusedrangeindex = 1
        return rs
    end

    if st == beforestartindex(rs)  # the index `i` is before the first range
        if rs.rstarts[1] - i > 1  # there is will be gap in indices after inserting
            pushfirst!(rs.rstarts, i)
            pushfirst!(rs.rstops, i)
        else  # prepend to first range
            rs.rstarts[1] = i
        end
        rs.lastusedrangeindex = 1
        return rs
    end

    r = get_range(rs, st)

    if i >= rs.rstarts[end]  # the index `i` is after the last range start
        if i > last(r) + 1  # there is will be the gap in indices after inserting
            push!(rs.rstarts, i)
            push!(rs.rstops, i)
        else  # just append to last range
            rs.rstops[st] += 1
        end
        rs.lastusedrangeindex = lastindex(rs)
        return rs
    end

    # the index `i` is somewhere between indices
    stnext = advance(rs, st)
    rnext = get_range(rs, stnext)

    if first(rnext) - last(r) == 2  # join ranges
        rs.rstarts[st] = first(r)
        rs.rstops[st] = last(rnext)
        deleteat!(rs.rstarts, stnext)
        deleteat!(rs.rstops, stnext)
        rs.lastusedrangeindex = st
    elseif i - last(r) == 1  # append to left range
        rs.rstops[st] += 1
        rs.lastusedrangeindex = st
    elseif first(rnext) - i == 1  # prepend to right range
        rs.rstarts[stnext] -= 1
        rs.lastusedrangeindex = stnext
    else  # insert single element chunk
        insert!(rs.rstarts, stnext, i)
        insert!(rs.rstops, stnext, i)
        rs.lastusedrangeindex = stnext
    end

    return rs

end


function Base.push!(rs::UnitRangesSortedSet{Ti}, idx) where {Ti}
    i = Ti(idx)

    st = check_exist_and_update(rs, i)
    st == nothing && return rs

    if length(rs) == 0
        rs.ranges[i] = i
        rs.lastusedrangeindex = firstindex(rs)
        return rs
    end

    if st == beforestartindex(rs)  # the index `i` is before the first index
        stnext = advance(rs, st)
        rnext = get_range(rs, stnext)
        if first(rnext) - i > 1  # there is will be gap in indices after inserting
            rs.ranges[i] = i
        else  # prepend to first range
            rs.ranges[i] = last(rnext)
            delete!((rs.ranges, stnext))
        end
        rs.lastusedrangeindex = firstindex(rs)
        return rs
    end

    r = get_range(rs, st)

    if i > last(last(rs)) # the index `i` is after the last range end index
        if last(r) + 1 < i  # there is will be the gap in indices after inserting
            rs.ranges[i] = i
        else  # just append to last range
            rs.ranges[st] = last(r)+1
        end
        rs.lastusedrangeindex = lastindex(rs)
        return rs
    end

    rs.lastusedrangeindex = beforestartindex(rs)

    # the index `i` is somewhere between indices
    stnext = advance(rs, st)
    rnext = get_range(rs, stnext)

    if first(rnext) - last(r) == 2  # join ranges
        rs.ranges[st] = last(rnext)
        delete!((rs.ranges, stnext))
    elseif i - last(r) == 1  # append to left range
        rs.ranges[st] = i
        rs.lastusedrangeindex = st
    elseif first(rnext) - i == 1  # prepend to right range
        rs.ranges[i] = last(rnext)
        delete!((rs.ranges, stnext))
    else  # insert single element range
        rs.ranges[i] = i
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
    todelete1 = get_range_start(rs, first(Iposition)) < first(I) ? first(Iposition) + 1 : first(Iposition)
    todelete2 = get_range_stop(rs, last(Iposition)) > last(I) ? last(Iposition) - 1 : last(Iposition)
    splice!(rs.rstarts, todelete1:todelete2)
    splice!(rs.rstops, todelete1:todelete2)

    # remains the side ranges to delete
    # Note: there is not possible `beforestartindex` or `pastendindex` in `Iposition`
    Iposition = searchsortedrange(rs, I)
    length(Iposition) == 0 && return rs

    Iposition_range_start = get_range_start(rs, first(Iposition))
    Iposition_range_stop = get_range_stop(rs, last(Iposition))

    # inside one range, thus split it
    if length(Iposition) == 1 &&
       Iposition_range_start < first(I) && last(I) < Iposition_range_stop
        insert!(rs.rstarts, advance(rs, first(Iposition)), last(I) + 1)
        insert!(rs.rstops, first(Iposition), first(I) - 1)

    # `I` is whole range, delete it
    elseif length(Iposition) == 1 &&
           Iposition_range_start >= first(I) && last(I) >= Iposition_range_stop
        deleteat!(rs.rstarts, first(Iposition))
        deleteat!(rs.rstops, first(Iposition))

    # `I` intersects with or inside in one range from left side, thus shrink from left side
    elseif length(Iposition) == 1 &&
           Iposition_range_start >= first(I) && last(I) < Iposition_range_stop
        rs.rstarts[first(Iposition)] = last(I) + 1

    # `I` intersects with or inside in one range from right side, thus shrink from right side
    # inside one range, ended to left side, thus shrink from right side
    elseif length(Iposition) == 1 &&
           Iposition_range_start < first(I) && last(I) >= Iposition_range_stop
        rs.rstops[first(Iposition)] = first(I) - 1

    # All remaining cases are with two ranges from `rs`
    elseif length(Iposition) != 2
        throw(AssertionError("FIXME: I Am an error in program logic."))

    else

        # delete whole second range
        if last(I) >= Iposition_range_stop
            deleteat!(rs.rstarts, last(Iposition))
            deleteat!(rs.rstops, last(Iposition))

        # or shrink it from left side
        else
            rs.rstarts[last(Iposition)] = last(I) + 1
        end

        # delete whole first range
        if Iposition_range_start >= first(I)
            deleteat!(rs.rstarts, first(Iposition))
            deleteat!(rs.rstops, first(Iposition))

        # or shrink it from right side
        else
            rs.rstops[first(Iposition)] = first(I) - 1
        end

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
    todelete1 = get_range_start(rs, first(Iposition)) < first(I) ? advance(rs, first(Iposition)) : first(Iposition)
    todelete2 = get_range_stop(rs, last(Iposition)) > last(I) ? regress(rs, last(Iposition)) : last(Iposition)
    for st in onlysemitokens(inclusive(rs.ranges, todelete1, todelete2))
        delete!((rs.ranges, st))
    end

    # remains the side ranges to delete
    # Note: there is not possible `beforestartindex` or `pastendindex` in `Iposition`
    Iposition = searchsortedrange(rs, I)
    length(Iposition) == 0 && return rs

    Iposition_range_start = get_range_start(rs, first(Iposition))
    Iposition_range_stop = get_range_stop(rs, last(Iposition))

    # inside one range, thus split it
    if length(Iposition) == 1 &&
       Iposition_range_start < first(I) && last(I) < Iposition_range_stop
        rs.ranges[last(I) + 1] = Iposition_range_stop
        rs.ranges[first(Iposition)] = first(I) - 1

    # `I` is whole range, delete it
    elseif length(Iposition) == 1 &&
           Iposition_range_start >= first(I) && last(I) >= Iposition_range_stop
        delete!((rs.ranges, first(Iposition)))

    # `I` intersects with or inside in one range from left side, thus shrink from left side
    elseif length(Iposition) == 1 &&
           Iposition_range_start >= first(I) && last(I) < Iposition_range_stop
        rs.ranges[last(I) + 1] = Iposition_range_stop
        delete!((rs.ranges, first(Iposition)))

    # `I` intersects with or inside in one range from right side, thus shrink from right side
    # inside one range, ended to left side, thus shrink from right side
    elseif length(Iposition) == 1 &&
           Iposition_range_start < first(I) && last(I) >= Iposition_range_stop
        rs.ranges[first(Iposition)] = first(I) - 1

    # All remaining cases are with two ranges from `rs`
    elseif length(Iposition) != 2
        throw(AssertionError("FIXME: I Am an error in program logic."))

    else

        # delete whole second range
        if last(I) >= Iposition_range_stop
            delete!((rs.ranges, last(Iposition)))

        # or shrink it from left side
        else
            rs.ranges[last(I) + 1] = Iposition_range_stop
            delete!((rs.ranges, last(Iposition)))
        end

        # delete whole first range
        if Iposition_range_start >= first(I)
            delete!((rs.ranges, first(Iposition)))

        # or shrink it from right side
        else
            rs.ranges[first(Iposition)] = first(I) - 1
        end

    end

    return rs
end


function Base.delete!(rs::UnitRangesSortedVector{Ti}, idx) where {Ti}
    i = Ti(idx)

    length(rs) == 0 && return rs

    st = searchsortedlastrange(rs, i)

    if st == beforestartindex(rs)  # the index `i` is before first range
        return rs
    end

    r = get_range(rs, st)

    if i > last(r)  # the index `i` is outside of range
        return rs
    end

    if last(r) - first(r) + 1 == 1
        deleteat!(rs.rstarts, st)
        deleteat!(rs.rstops, st)
    elseif i == last(r)  # last index in range
        rs.rstops[st] -= 1
    elseif i == first(r)  # first index in range
        rs.rstarts[st] += 1
    else
        insert!(rs.rstarts, advance(rs, st), i+1)
        insert!(rs.rstops, advance(rs, st), last(r))
        rs.rstops[st] = i-1
    end

    rs.lastusedrangeindex = beforestartindex(rs)

    return rs
end


function Base.delete!(rs::UnitRangesSortedSet{Ti}, idx) where {Ti}
    i = Ti(idx)

    length(rs) == 0 && return rs

    st = searchsortedlastrange(rs, i)

    if st == beforestartindex(rs)  # the index `i` is before first index
        return rs
    end

    r = get_range(rs, st)

    if i > last(r)  # the index `i` is outside of range
        return rs
    end

    if last(r) - first(r) + 1 == 1
        delete!((rs.ranges, st))
    elseif i == last(r)  # last index in range
        rs.ranges[st] = last(r) - 1
    elseif i == first(r)  # first index in range
        delete!((rs.ranges, st))
        rs.ranges[i+1] = last(r)
    else
        rs.ranges[st] = i - 1
        rs.ranges[i+1] = last(r)
    end

    rs.lastusedrangeindex = beforestartindex(rs)

    return rs
end

function Base.pop!(rs::AbstractUnitRangesSortedSet, i)
    if in(i, rs)
        delete!(rs, i)
        return i
    else
        throw(KeyError(i))
    end
end
function Base.pop!(rs::AbstractUnitRangesSortedSet, i, default)
    if in(i, rs)
        delete!(rs, i)
        return i
    else
        return default
    end
end


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

# TODO: intersect (and view(rs, r::InitRange)), setdiff, symdiff, issubset, eachindex (SomeEachindexIterator), map
# setdiff already exist in julia/base/abstractset.jl

Base.copy(rs::T) where {T<:SubUnitRangesSortedSet} = T(rs.parent, rs.start, rs.stop)
Base.copy(rs::T) where {T<:UnitRangesSortedVector} = T(rs.lastusedrangeindex, copy(rs.rstarts), copy(rs.rstops))
Base.copy(rs::T) where {T<:UnitRangesSortedSet} = T(rs.lastusedrangeindex, packcopy(rs.ranges))
Base.copymutable(rs::AbstractUnitRangesSortedSet) = copy(rs)

Base.emptymutable(rs::T) where {T<:AbstractUnitRangesSortedSet} = T()


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
        issubset(r1, rs2) || return false
    end
    return true
end
function Base.issubset(rs1::Union{AbstractSet,AbstractVector,AbstractRange,Tuple}, rs2::AbstractUnitRangesSortedSet)
    for r1 in rs1
        issubset(r1, rs2) || return false
    end
    return true
end
function Base.issubset(rs1::AbstractUnitRangesSortedSet, rs2::Union{AbstractSet,AbstractVector,AbstractRange,Tuple})
    for r1 in rs1
        issubset(r1, rs2) || return false
    end
    return true
end


function Base.filter(pred::Function, rs::AbstractUnitRangesSortedSet)
    res = Base.emptymutable(rs)
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
#
#  Aux functions
#

Base.show(io::IO, ur::URSSUnitRange) = print(io, repr(first(ur)), ':', length(ur) == 0 ? "0" : "1" , ':', repr(last(ur)))

function Base.show(io::IO, ::MIME"text/plain", x::URSSUnitRange)
    len = length(x)
    print(io, len, "-element range ", first(x), ":", last(x), " for ", typeof(x.parent))
    if len != 0
        println(io, " containing:")
        show(IOContext(io, :typeinfo => eltype(x)), x)
    end
end

function Base.show(io::IOContext, x::URSSUnitRange)
    ranges = [get_range(x, s) for s in x]

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
