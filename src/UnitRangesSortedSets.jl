
module UnitRangesSortedSets
export AbstractUnitRangesSortedSet, UnitRangesSortedVector, UnitRangesSortedSet, SubUnitRangesSortedSet
export findfirstnz, findlastnz, findfirstnzindex, findlastnzindex
export iteratenzpairs, iteratenzpairsview, iteratenzvalues, iteratenzvaluesview, iteratenzindices
export testfun_create, testfun_createSV, testfun_createVL, testfun_create_seq, testfun_create_dense, testfun_delete!, testfun_in, testfun_in_outer, testfun_in_rand, testfun_in_seq, testfun_nzgetindex, testfun_setindex!, testfun_nzchunks, testfun_nzpairs, testfun_nzindices, testfun_nzvalues, testfun_nzvaluesview, testfun_findnz
export searchsortedrange, searchsortedfirstrange, searchsortedlastrange


import Base: ForwardOrdering, Forward
const FOrd = ForwardOrdering

using DocStringExtensions
using DataStructures
using FillArrays
using IterTools
using Setfield
using SparseArrays
using StaticArrays
import SparseArrays: indtype, nonzeroinds, nonzeros
using Random


abstract type AbstractUnitRangesSortedSet{Ti} <: AbstractSet{Ti} end

struct SubUnitRangesSortedSet{Ti,P} <: AbstractUnitRangesSortedSet{Ti}
    parent::P
    start::Ti
    stop::Ti
    startindex::Int
    stopindex::Int
end

Base.@propagate_inbounds function Base.view(rs::Tp, I::UnitRange) where {Tp<:AbstractUnitRangesSortedSet{Ti}} where Ti
    r = searchsortedrange(rs, I)
    return SubUnitRangesSortedSet{Ti,Tp}(rs, Ti(I.start), Ti(I.stop), r.start, r.stop)
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

function UnitRangesSortedVector(values::Union{AbstractVector, AbstractSet, Tuple})
    rs = UnitRangesSortedVector{eltype(values)}()
    for r in values
        push!(rs, r)
    end
    rs
end
function UnitRangesSortedVector{Ti}(values::Union{AbstractVector, AbstractSet, Tuple}) where {Ti}
    rs = UnitRangesSortedVector{Ti}()
    for r in values
        push!(rs, r)
    end
    rs
end

function UnitRangesSortedVector(rs::AbstractUnitRangesSortedSet{Ti}) where {Ti}
    rstarts = Vector{Ti}(undef, length(rs))
    rstops = Vector{Ti}(undef, length(rs))
    for (i, r) in enumerate(rs)
        rstarts[i] = r.start
        rstops[i] = r.stop
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

function UnitRangesSortedSet(values::Union{AbstractVector, AbstractSet, Tuple})
    rs = UnitRangesSortedSet{eltype(values)}()
    for r in values
        push!(rs, r)
    end
    rs
end
function UnitRangesSortedSet{Ti}(values::Union{AbstractVector, AbstractSet, Tuple}) where {Ti}
    rs = UnitRangesSortedSet{Ti}()
    for r in values
        push!(rs, r)
    end
    rs
end

function UnitRangesSortedSet(rs::AbstractUnitRangesSortedSet{Ti}) where {Ti}
    ranges = SortedDict{Ti,Ti,FOrd}(Forward)
    for r in rs
        ranges[r.start] = r.stop
    end
    UnitRangesSortedSet{Ti}(beforestartsemitoken(ranges), ranges)
end




Base.length(rs::UnitRangesSortedVector) = length(rs.rstarts)
Base.length(rs::UnitRangesSortedSet) = length(rs.ranges)
Base.haslength(rs::AbstractUnitRangesSortedSet) = true
Base.hasfastin(rs::AbstractUnitRangesSortedSet) = true
Base.isempty(rs::AbstractUnitRangesSortedSet) = length(rs) == 0
Base.size(rs::AbstractUnitRangesSortedSet) = (length(rs),)
#Base.axes(rs::AbstractUnitRangesSortedSet) = (Base.OneTo(rs.n),)
Base.eltype(::AbstractUnitRangesSortedSet{Ti}) where {Ti} = Ti
#Base.IndexStyle(::AbstractUnitRangesSortedSet) = IndexLinear()


#function Base.collect(::Type{ElType}, rs::AbstractUnitRangesSortedSet) where ElType
#    res = Vector{ElType}(undef, foldl((s,r)->s+(r.stop-r.start+1), rs, init=0))
#    i = 0
#    for r in rs, el in r
#        res[i+=1] = ElType(el)
#    end
#    return res
#end
#Base.collect(rs::AbstractUnitRangesSortedSet) = collect(eltype(rs), rs)
function Base.collect(::Type{ElType}, rs::AbstractUnitRangesSortedSet) where ElType
    res = Vector{UnitRange{ElType}}(undef, length(rs))
    i = 0
    for r in rs
        res[i+=1] = ElType(r.start):ElType(r.stop)
    end
    return res
end
function Base.collect(rs::AbstractUnitRangesSortedSet{Ti}) where Ti
    res = Vector{UnitRange{Ti}}(undef, length(rs))
    i = 0
    for r in rs
        res[i+=1] = r.start:r.stop
    end
    return res
end

Base.@propagate_inbounds length_of_that_range(rs::AbstractUnitRangesSortedSet, r) = r.stop - r.start + 1
@inline get_range_length(rs::AbstractUnitRangesSortedSet, i) = ((start, stop) = get_range_indices(rs, i); stop - start + 1)
@inline get_range_indices(rs::UnitRangesSortedVector, i) = (rs.rstarts[i], rs.rstops[i])
@inline get_range_indices(rs::UnitRangesSortedSet, i::DataStructures.Tokens.IntSemiToken) = deref((rs.ranges, i))
@inline get_range(rs::UnitRangesSortedVector, i) = (rs.rstarts[i]:rs.rstops[i])
@inline get_range(rs::UnitRangesSortedSet, i::DataStructures.Tokens.IntSemiToken) = ((start, stop) = deref((rs.ranges, i)); return (start:stop))
@inline function get_range(rs::SubArray{<:Any,<:Any,<:T}, i) where {T<:AbstractUnitRangesSortedSet}
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
@inline function get_range(rs::SubArray{<:Any,<:Any,<:T}, i) where {T<:AbstractUnitRangesSortedSet}
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
@inline Base.first(rs::UnitRangesSortedVector) = (rs.rstarts[1]:rs.rstops[1])
@inline Base.first(rs::UnitRangesSortedSet) = ((start,stop) = deref((rs.ranges, startof(rs.ranges))); (start:stop))
@inline Base.last(rs::UnitRangesSortedVector) = (rs.rstarts[end]:rs.rstops[end])
@inline Base.last(rs::UnitRangesSortedSet) = ((start,stop) = deref((rs.ranges, lastindex(rs.ranges))); (start:stop))

@inline beforestartindex(rs::UnitRangesSortedVector) = firstindex(rs.rstarts) - 1
@inline beforestartindex(rs::UnitRangesSortedSet) = beforestartsemitoken(rs.ranges)
@inline pastendindex(rs::UnitRangesSortedVector) = lastindex(rs.rstarts) + 1
@inline pastendindex(rs::UnitRangesSortedSet) = pastendsemitoken(rs.ranges)

@inline DataStructures.advance(rs::UnitRangesSortedVector, state) = state + 1
@inline DataStructures.advance(rs::UnitRangesSortedSet, state) = advance((rs.ranges, state))

@inline urangevectorisless(x::UnitRange, y::UnitRange) = x.start < y.start
@inline urangevectorisless(x::UnitRange, y) = x.start < y
@inline urangevectorisless(x, y::UnitRange) = x < y.start

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
@inline searchsortedrange(rs::AbstractUnitRangesSortedSet, I::UnitRange) =
    searchsortedfirstrange(rs, I.start):searchsortedlastrange(rs, I.stop)
#@inline searchsortedrange(rs::UnitRangesSortedSet, I::UnitRange) =
#    (searchsortedfirstrange(rs, I.start), searchsortedlastrange(rs, I.stop))

"Returns indexes of range in `rs` in which `i` may be inserted. Or negative range in the case of `i` is
 between `rs` ranges, and indices of resulted range is the indexes of that neighbors."
@inline function searchsortedrange(rs::AbstractUnitRangesSortedSet, i)
    st = searchsortedlastrange(rs, i)
    if st != beforestartindex(rs) && in(i, get_range(rs, st))
        return st:st
    else
        return advance(rs, st):st
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

Base.@propagate_inbounds index_status(rs::UnitRangesSortedVector, st) = 0 <= st <= length(rs)+1 ? 1 : 0
Base.@propagate_inbounds index_status(rs::UnitRangesSortedSet, st) = status((rs.ranges, st))

"Returns the index of first non-zero element in sparse vector."
@inline findfirstnzindex(rs::SparseVector) = nnz(rs) > 0 ? rs.ranges[1] : nothing
@inline findfirstnzindex(rs::AbstractUnitRangesSortedSet{Ti}) where {Ti} =
    nnz(rs) > 0 ? Ti(rs.ranges[1]) : nothing
@inline findfirstnzindex(rs::AbstractUnitRangesSortedSet{Ti}) where {Ti} =
    nnz(rs) > 0 ? Ti(deref_key((rs.ranges, startof(rs.ranges)))) : nothing
function findfirstnzindex(rs::SubArray{<:Any,<:Any,<:T})  where {T<:AbstractUnitRangesSortedSet{Ti}} where {Ti}
    nnz(rs.parent) == 0 && return nothing
    ifirst, r_stop = first(rs.indices[1]), last(rs.indices[1])
    st = searchsortedlast_nzchunk(rs.parent, ifirst)
    st == pastendindex(rs.parent) && return nothing
    key = get_range_start(rs.parent, st)
    len = get_range_length(rs.parent, st)
    if key <= ifirst < key + len  # ifirst index within nzchunk range
        return Ti(1)
    elseif ifirst <= key <= r_stop  # nzchunk[1] somewhere in ifirst:r_stop range
        return Ti(key-ifirst+1)
    else
        return nothing
    end
end

"Returns the index of last non-zero element in sparse vector."
@inline findlastnzindex(rs::SparseVector) = nnz(rs) > 0 ? rs.ranges[end] : nothing
@inline findlastnzindex(rs::AbstractUnitRangesSortedSet) =
    nnz(rs) > 0 ? rs.ranges[end] + length_of_that_range(rs, rs.ranges[end]) - 1 : nothing
@inline function findlastnzindex(rs::AbstractUnitRangesSortedSet)
    if nnz(rs) > 0
        lasttoken = lastindex(rs.ranges)
        return deref_key((rs.ranges, lasttoken)) + length_of_that_range(rs, deref_value((rs.ranges, lasttoken))) - 1
    else
        return nothing
    end
end
function findlastnzindex(rs::SubArray{<:Any,<:Any,<:T})  where {T<:AbstractUnitRangesSortedSet{Ti}} where {Ti}
    nnz(rs.parent) == 0 && return nothing
    ifirst, r_stop = first(rs.indices[1]), last(rs.indices[1])
    st = searchsortedfirst_nzchunk(rs.parent, r_stop)
    st == beforestartindex(rs.parent) && return nothing
    key = get_range_start(rs.parent, st)
    len = get_range_length(rs.parent, st)
    if key <= r_stop < key + len  # r_stop index within nzchunk range
        return Ti(r_stop - ifirst + 1)
    elseif ifirst <= key+len-1 <= r_stop  # nzchunk[end] somewhere in ifirst:r_stop range
        return Ti(key+len-1 - ifirst+1)
    else
        return nothing
    end
end

"Returns value of first non-zero element in the sparse vector."
@inline findfirstnz(rs::AbstractSparseVector) = nnz(rs) > 0 ? rs[findfirstnzindex(rs)] : nothing
function findfirstnz(rs::SubArray{<:Any,<:Any,<:T})  where {T<:AbstractUnitRangesSortedSet{Ti}} where {Ti}
    nnz(rs.parent) == 0 && return nothing
    ifirst, r_stop = first(rs.indices[1]), last(rs.indices[1])
    st = searchsortedlast_nzchunk(rs.parent, ifirst)
    st == pastendindex(rs.parent) && return nothing
    key, chunk = get_range(rs.parent, st)
    len = length_of_that_range(rs.parent, chunk)
    if key <= ifirst < key + len  # ifirst index within nzchunk range
        return chunk[ifirst-key+1]
    elseif ifirst <= key <= r_stop  # nzchunk[1] somewhere in ifirst:r_stop range
        return chunk[1]
    else
        return nothing
    end
end

"Returns value of last non-zero element in the sparse vector."
@inline findlastnz(rs::AbstractSparseVector) = nnz(rs) > 0 ? rs[findlastnzindex(rs)] : nothing
function findlastnz(rs::SubArray{<:Any,<:Any,<:T})  where {T<:AbstractUnitRangesSortedSet{Ti}} where {Ti}
    nnz(rs.parent) == 0 && return nothing
    ifirst, r_stop = first(rs.indices[1]), last(rs.indices[1])
    st = searchsortedfirst_nzchunk(rs.parent, r_stop)
    st == beforestartindex(rs.parent) && return nothing
    key, chunk = get_range(rs.parent, st)
    len = length_of_that_range(rs.parent, chunk)
    if key <= r_stop < key + len  # r_stop index within nzchunk range
        return chunk[r_stop-key+1]
    elseif ifirst <= key+len-1 <= r_stop  # nzchunk[end] somewhere in ifirst:r_stop range
        return chunk[end]
    else
        return nothing
    end
end


@inline function Base.findfirst(testf::Function, rs::AbstractUnitRangesSortedSet)
    for p in nzpairs(rs)
        testf(last(p)) && return first(p)
    end
    return nothing
end

@inline Base.findall(testf::Function, rs::AbstractUnitRangesSortedSet) = collect(first(p) for p in nzpairs(rs) if testf(last(p)))


#
#  Iterators
#


Base.@propagate_inbounds function Base.iterate(rs::AbstractUnitRangesSortedSet, state = firstindex(rs))
    if state != pastendindex(rs)
        return (get_range(rs, state), advance(rs, state))
    else
        return nothing
    end
end

#
# Assignments
#


#@inline Base.in(II::UnitRange, rs::AbstractUnitRangesSortedSet{Ti}) where Ti = in(convert(UnitRange{Ti}, II), rs)

@inline function Base.in(II::UnitRange{Ti}, rs::AbstractUnitRangesSortedSet{Ti}) where Ti
    I = convert(UnitRange{Ti}, II)
    # fast check for cached range index
    if (st = rs.lastusedrangeindex) != beforestartindex(rs)
        if issubset(I, get_range(rs, st))
            return true
        end
    end
    # cached range index miss (or index not stored), thus try search
    Irange = searchsortedrange(rs, I)
    rs.lastusedrangeindex = Irange.stop
    if length(Irange) != 1
        return false
    elseif issubset(I, get_range(rs, Irange.start))
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


@inline Base.haskey(rs::AbstractUnitRangesSortedSet, idx) = in(idx, rs)

function Base.push!(rs::UnitRangesSortedVector{Ti}, II::UnitRange) where Ti
    I = convert(UnitRange{Ti}, II)

    if length(I) < 2
        for i in I
            push!(rs, i)
        end
        return rs
    end

    rs.lastusedrangeindex = beforestartindex(rs)

    Irange = searchsortedrange(rs, I)
    if length(Irange) == 1 && I == get_range(rs, Irange.start)
        # `I` already exist in `rs` and the same as the corresponding range in `rs`, nothing to do
        return rs
    else
        # delete all inside `I` range in `rs`
        delete!(rs, I)
    end

    # get neighbors pointers; range will be zero-length
    # Note: `Irange.start` will be point to range in `rs` on right side to `I` and
    # `Irange.stop` will be point to range in `rs` on left side to `I`.
    Irange = searchsortedrange(rs, I.start)
    @boundscheck @assert Irange == searchsortedrange(rs, I.stop) "FIXME: I Am an error in program logic."

    # `I` is adjoined in both sides with `rs` ranges, thus join them all
    if Irange.stop != beforestartindex(rs) && get_range_stop(rs, Irange.stop) + 1 == I.start &&
       Irange.start != pastendindex(rs) && get_range_start(rs, Irange.start) - 1 == I.stop
        rs.rstops[Irange.stop] = get_range_stop(rs, Irange.start)
        deleteat!(rs.rstarts, Irange.start)
        deleteat!(rs.rstops, Irange.start)

    # `I` is adjoin with `rs` range on left side, thus append to left range
    elseif Irange.stop != beforestartindex(rs) && get_range_stop(rs, Irange.stop) + 1 == I.start &&
           Irange.start != pastendindex(rs) && get_range_start(rs, Irange.start) - 1 != I.stop
        rs.rstops[Irange.stop] = I.stop

    # `I` is adjoin with `rs` range on right side, thus prepend to right range
    elseif Irange.stop != beforestartindex(rs) && get_range_stop(rs, Irange.stop) + 1 != I.start &&
           Irange.start != pastendindex(rs) && get_range_start(rs, Irange.start) - 1 == I.stop
        rs.rstarts[Irange.start] = I.start

    # `I` is separate from both sides, insert it
    elseif Irange.stop != beforestartindex(rs) && get_range_stop(rs, Irange.stop) + 1 != I.start &&
           Irange.start != pastendindex(rs) && get_range_start(rs, Irange.start) - 1 != I.stop
        insert!(rs.rstarts, Irange.start, I.start)
        insert!(rs.rstops, Irange.start, I.stop)

    # `I` is first range in `rs` and adjoin with range on right side, thus prepend `I` to right range
    elseif Irange.stop == beforestartindex(rs) &&
           Irange.start != pastendindex(rs) && get_range_start(rs, Irange.start) - 1 == I.stop
        rs.rstarts[Irange.start] = I.start

    # `I` is separate first range in `rs`, insert it
    elseif Irange.stop == beforestartindex(rs) &&
           Irange.start != pastendindex(rs) && get_range_start(rs, Irange.start) - 1 != I.stop
        insert!(rs.rstarts, Irange.start, I.start)
        insert!(rs.rstops, Irange.start, I.stop)

    # `I` is last range in `rs` and adjoin with range on left side, thus append `I` to left range
    elseif Irange.stop != beforestartindex(rs) && get_range_stop(rs, Irange.stop) + 1 == I.start &&
           Irange.start == pastendindex(rs)
        rs.rstops[Irange.stop] = I.stop

    # `I` is separate last range in `rs`, insert it
    elseif Irange.stop != beforestartindex(rs) && get_range_stop(rs, Irange.stop) + 1 != I.start &&
           Irange.start == pastendindex(rs)
        insert!(rs.rstarts, Irange.start, I.start)
        insert!(rs.rstops, Irange.start, I.stop)

    # `rs` is empty
    else #if Irange.stop == beforestartindex(rs) && Irange.start == pastendindex(rs)
        insert!(rs.rstarts, 1, I.start)
        insert!(rs.rstops, 1, I.stop)
    end

    return rs
end

function Base.push!(rs::UnitRangesSortedSet{Ti}, II::UnitRange) where Ti
    I = convert(UnitRange{Ti}, II)

    if length(I) < 2
        for i in I
            push!(rs, i)
        end
        return rs
    end

    rs.lastusedrangeindex = beforestartindex(rs)

    Irange = searchsortedrange(rs, I)
    if length(Irange) == 1 && I == get_range(rs, Irange.start)
        # `I` already exist in `rs` and the same as the corresponding range in `rs`, nothing to do
        return rs
    else
        # delete all inside `I` range in `rs`
        delete!(rs, I)
    end

    # get neighbors pointers; range will be zero-length
    # Note: `Irange.start` will be point to range in `rs` on right side to `I` and
    # `Irange.stop` will be point to range in `rs` on left side to `I`.
    Irange = searchsortedrange(rs, I.start)
    @boundscheck @assert Irange == searchsortedrange(rs, I.stop) "FIXME: I Am an error in program logic."

    # `I` is adjoined in both sides with `rs` ranges, thus join them all
    if Irange.stop != beforestartindex(rs) && get_range_stop(rs, Irange.stop) + 1 == I.start &&
       Irange.start != pastendindex(rs) && get_range_start(rs, Irange.start) - 1 == I.stop
        rs.rstops[Irange.stop] = get_range_stop(rs, Irange.start)
        deleteat!(rs.rstarts, Irange.start)
        deleteat!(rs.rstops, Irange.start)

    # `I` is adjoin with `rs` range on left side, thus append to left range
    elseif Irange.stop != beforestartindex(rs) && get_range_stop(rs, Irange.stop) + 1 == I.start &&
           Irange.start != pastendindex(rs) && get_range_start(rs, Irange.start) - 1 != I.stop
        rs.rstops[Irange.stop] = I.stop

    # `I` is adjoin with `rs` range on right side, thus prepend to right range
    elseif Irange.stop != beforestartindex(rs) && get_range_stop(rs, Irange.stop) + 1 != I.start &&
           Irange.start != pastendindex(rs) && get_range_start(rs, Irange.start) - 1 == I.stop
        rs.rstarts[Irange.start] = I.start

    # `I` is separate from both sides, insert it
    elseif Irange.stop != beforestartindex(rs) && get_range_stop(rs, Irange.stop) + 1 != I.start &&
           Irange.start != pastendindex(rs) && get_range_start(rs, Irange.start) - 1 != I.stop
        insert!(rs.rstarts, Irange.start, I.start)
        insert!(rs.rstops, Irange.start, I.stop)

    # `I` is first range in `rs` and adjoin with range on right side, thus prepend `I` to right range
    elseif Irange.stop == beforestartindex(rs) &&
           Irange.start != pastendindex(rs) && get_range_start(rs, Irange.start) - 1 == I.stop
        rs.rstarts[Irange.start] = I.start

    # `I` is separate first range in `rs`, insert it
    elseif Irange.stop == beforestartindex(rs) &&
           Irange.start != pastendindex(rs) && get_range_start(rs, Irange.start) - 1 != I.stop
        insert!(rs.rstarts, Irange.start, I.start)
        insert!(rs.rstops, Irange.start, I.stop)

    # `I` is last range in `rs` and adjoin with range on left side, thus append `I` to left range
    elseif Irange.stop != beforestartindex(rs) && get_range_stop(rs, Irange.stop) + 1 == I.start &&
           Irange.start == pastendindex(rs)
        rs.rstops[Irange.stop] = I.stop

    # `I` is separate last range in `rs`, insert it
    elseif Irange.stop != beforestartindex(rs) && get_range_stop(rs, Irange.stop) + 1 != I.start &&
           Irange.start == pastendindex(rs)
        insert!(rs.rstarts, Irange.start, I.start)
        insert!(rs.rstops, Irange.start, I.stop)

    # `rs` is empty
    else #if Irange.stop == beforestartindex(rs) && Irange.start == pastendindex(rs)
        insert!(rs.rstarts, 1, I.start)
        insert!(rs.rstops, 1, I.stop)
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
    @boundscheck if index_status(rs, st) == 0 # invalid semitoken
        throw(KeyError(i))
    end

    # check the index exist and update its data
    if st != beforestartindex(rs)  # the index `i` is not before the first index
        if i <= last(get_range_indices(rs, st))
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
        if i > r.stop + 1  # there is will be the gap in indices after inserting
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

    if rnext.start - r.stop == 2  # join ranges
        rs.rstarts[st] = r.start
        rs.rstops[st] = rnext.stop
        deleteat!(rs.rstarts, stnext)
        deleteat!(rs.rstops, stnext)
        rs.lastusedrangeindex = st
    elseif i - r.stop == 1  # append to left range
        rs.rstops[st] += 1
        rs.lastusedrangeindex = st
    elseif rnext.start - i == 1  # prepend to right range
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
        if rnext.start - i > 1  # there is will be gap in indices after inserting
            rs.ranges[i] = i
        else  # prepend to first range
            rs.ranges[i] = rnext.stop
            delete!((rs.ranges, stnext))
        end
        rs.lastusedrangeindex = firstindex(rs)
        return rs
    end

    r = get_range(rs, st)

    if i > last(rs).stop # the index `i` is after the last range end index
        if r.stop + 1 < i  # there is will be the gap in indices after inserting
            rs.ranges[i] = i
        else  # just append to last range
            rs.ranges[st] = r.stop+1
        end
        rs.lastusedrangeindex = lastindex(rs)
        return rs
    end

    rs.lastusedrangeindex = beforestartindex(rs)

    # the index `i` is somewhere between indices
    stnext = advance(rs, st)
    rnext = get_range(rs, stnext)

    if rnext.start - r.stop == 2  # join ranges
        rs.ranges[st] = rnext.stop
        delete!((rs.ranges, stnext))
    elseif i - r.stop == 1  # append to left range
        rs.ranges[st] = i
        rs.lastusedrangeindex = st
    elseif rnext.start - i == 1  # prepend to right range
        rs.ranges[i] = rnext.stop
        delete!((rs.ranges, stnext))
    else  # insert single element range
        rs.ranges[i] = i
    end

    return rs

end


function Base.delete!(rs::UnitRangesSortedVector{Ti}, II::UnitRange) where {Ti}
    length(II) == 0 && return rs

    I = convert(UnitRange{Ti}, II)

    Irange = searchsortedrange(rs, I)
    length(Irange) == 0 && return rs

    rs.lastusedrangeindex = beforestartindex(rs)

    # delete all concluded ranges
    todelete1 = get_range_start(rs, Irange.start) < I.start ? Irange.start + 1 : Irange.start
    todelete2 = get_range_stop(rs, Irange.stop) > I.stop ? Irange.stop - 1 : Irange.stop
    splice!(rs.rstarts, todelete1:todelete2)
    splice!(rs.rstops, todelete1:todelete2)

    # remains the side ranges to delete
    # Note: there is not possible `beforestartindex` or `pastendindex` in `Irange`
    Irange = searchsortedrange(rs, I)
    length(Irange) == 0 && return rs

    # inside one range, thus split it
    if length(Irange) == 1 &&
       get_range_start(rs, Irange.start) < I.start && I.stop < get_range_stop(rs, Irange.stop)
        insert!(rs.rstarts, advance(rs, Irange.start), I.stop + 1)
        insert!(rs.rstops, Irange.start, I.start - 1)

    # `I` is whole range, delete it
    elseif length(Irange) == 1 &&
           get_range_start(rs, Irange.start) >= I.start && I.stop >= get_range_stop(rs, Irange.stop)
        deleteat!(rs.rstarts, Irange.start)
        deleteat!(rs.rstops, Irange.start)

    # `I` intersects with or inside in one range from left side, thus shrink from left side
    elseif length(Irange) == 1 &&
           get_range_start(rs, Irange.start) >= I.start && I.stop < get_range_stop(rs, Irange.stop)
        rs.rstarts[Irange.start] = I.stop + 1

    # `I` intersects with or inside in one range from right side, thus shrink from right side
    # inside one range, ended to left side, thus shrink from right side
    elseif length(Irange) == 1 &&
           get_range_start(rs, Irange.start) < I.start && I.stop >= get_range_stop(rs, Irange.stop)
        rs.rstops[Irange.start] = I.start - 1

    # All remaining cases are with two ranges from `rs`
    elseif length(Irange) != 2
        throw(AssertionError("FIXME: I Am an error in program logic."))

    else

        # delete whole second range
        if I.stop >= get_range_stop(rs, Irange.stop)
            deleteat!(rs.rstarts, Irange.stop)
            deleteat!(rs.rstops, Irange.stop)
        # or shrink it from left side
        else
            rs.rstarts[Irange.stop] = I.stop + 1
        end

        # delete whole first range
        if get_range_start(rs, Irange.start) >= I.start
            deleteat!(rs.rstarts, Irange.start)
            deleteat!(rs.rstops, Irange.start)
        # or shrink it from right side
        else
            rs.rstops[Irange.start] = I.start - 1
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

    if i > r.stop  # the index `i` is outside of range
        return rs
    end

    if r.stop - r.start + 1 == 1
        deleteat!(rs.rstarts, st)
        deleteat!(rs.rstops, st)
    elseif i == r.stop  # last index in range
        rs.rstops[st] -= 1
    elseif i == r.start  # first index in range
        rs.rstarts[st] += 1
    else
        insert!(rs.rstarts, advance(rs, st), i+1)
        insert!(rs.rstops, advance(rs, st), r.stop)
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

    if i > r.stop  # the index `i` is outside of range
        return rs
    end

    if r.stop - r.start + 1 == 1
        delete!((rs.ranges, st))
    elseif i == r.stop  # last index in range
        rs.ranges[st] = r.stop - 1
    elseif i == r.start  # first index in range
        delete!((rs.ranges, st))
        rs.ranges[i+1] = r.stop
    else
        rs.ranges[st] = i - 1
        rs.ranges[i+1] = r.stop
    end

    rs.lastusedrangeindex = beforestartindex(rs)

    return rs
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
Base.copy(rs::T) where {T<:UnitRangesSortedSet} = T(rs.lastusedrangeindex, copy(rs.ranges))
Base.copymutable(rs::AbstractUnitRangesSortedSet) = copy(rs)

Base.emptymutable(rs::T) where {T<:AbstractUnitRangesSortedSet} = T()


Base.union(rs::AbstractUnitRangesSortedSet, rss...) = union!(copy(rs), rss...)
Base.union(rs::AbstractUnitRangesSortedSet, rs2) = union!(copy(rs), rs2)

@inline Base.union!(rs::AbstractUnitRangesSortedSet, rss...) = union!(union!(rs, rss[1]), Iterators.tail(rss)...)
function Base.union!(rs::AbstractUnitRangesSortedSet, rs2)
    issubset(rs2, rs) && return rs
    for r in rs2
        push!(rs, r)
    end
    return rs
end

@inline Base.union!(rs::AbstractUnitRangesSortedSet, r::UnitRange) = push!(rs, r)

function Base.isequal(rs1::AbstractUnitRangesSortedSet, rs2::AbstractUnitRangesSortedSet)
    length(rs1) == length(rs2) || return false
    for (r1,r2) in zip(rs1, rs2)
        r1 == r2 || return false
    end
    return true
end

#function Base.issubset(rs1::AbstractUnitRangesSortedSet, rs2::AbstractUnitRangesSortedSet)
#    length(rs1) == length(rs2) || return false
#    for (r1,r2) in zip(rs1, rs2)
#        r1 == r2 || return false
#    end
#    return true
#end


##
##  Aux functions
##
## derived from stdlib/SparseArrays/src/sparsevector.jl
##
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
                print(io, lpad(repr(ranges[k].start), pad), ":", repr(ranges[k].stop))
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
            pad = max(pad, length(repr(r.start)), length(repr(r.stop)))
            #pad = max(pad, ndigits(r.start), ndigits(r.stop))
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


function testfun_nzchunks(sv)
    I = 0
    S = 0.0
    for (startindex,chunk) in nzchunkspairs(sv)
        startindex -= 1
        for i in axes(chunk,1)
            I += startindex + i
            S += chunk[i]
        end
    end
    (I, S)
end

function testfun_nzpairs(sv)
    I = 0
    S = 0.0
    for (k,v) in nzpairs(sv)
        I += k
        S += v
    end
    (I, S)
end

function testfun_nzindices(sv)
    I = 0
    for k in nzindices(sv)
        I += k
    end
    (I, 0.0)
end

function testfun_nzvalues(sv)
    S = 0.0
    for v in nzvalues(sv)
        S += v
    end
    (0, S)
end

function testfun_nzvaluesview(sv)
    S = 0.0
    for v in nzvaluesview(sv)
        S += v[1]
    end
    (0, S)
end

function testfun_findnz(sv)
    I = 0
    S = 0.0
    for (k,v) in zip(SparseArrays.findnz(sv)...)
        I += k
        S += v
    end
    (I, S)
end


end  # of module DensedSparseVectors
