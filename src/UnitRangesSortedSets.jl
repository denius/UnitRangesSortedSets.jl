
module UnitRangesSortedSets
export AbstractUnitRangesSortedSet, UnitRangesSortedVector, UnitRangesSortedSet, SubUnitRangesSortedSet
export findfirstnz, findlastnz, findfirstnzindex, findlastnzindex
export iteratenzpairs, iteratenzpairsview, iteratenzvalues, iteratenzvaluesview, iteratenzindices
export testfun_create, testfun_createSV, testfun_createVL, testfun_create_seq, testfun_create_dense, testfun_delete!, testfun_in, testfun_in_outer, testfun_in_rand, testfun_in_seq, testfun_nzgetindex, testfun_setindex!, testfun_nzchunks, testfun_nzpairs, testfun_nzindices, testfun_nzvalues, testfun_nzvaluesview, testfun_findnz


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
end

Base.@propagate_inbounds function Base.view(rs::Tp, I::UnitRange) where {Tp<:AbstractUnitRangesSortedSet{Ti}} where Ti
    return SubUnitRangesSortedSet{Ti,Tp}(rs, Ti(I.start), Ti(I.stop))
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
    ranges::Vector{UnitRange{Ti}}
end

UnitRangesSortedVector{Ti}() where {Ti} = UnitRangesSortedVector{Ti}(0, Vector{UnitRange{Ti}}(undef, 0))

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
    ranges = Vector{UnitRange{Ti}}(undef, length(rs))
    for (i, r) in enumerate(rs)
        ranges[i] = r
    end
    UnitRangesSortedVector{Ti}(firstindex(ranges) - 1, ranges)
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




Base.length(rs::AbstractUnitRangesSortedSet) = length(rs.ranges)
Base.isempty(rs::AbstractUnitRangesSortedSet) = length(rs) == 0
#Base.size(rs::AbstractUnitRangesSortedSet) = (rs.n,)
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

#Base.@propagate_inbounds length_of_that_range(rs::UnitRangesSortedVector, chunk) = chunk
#Base.@propagate_inbounds length_of_that_range(rs::UnitRangesSortedSet, chunk) = chunk
@inline get_range_length(rs::AbstractUnitRangesSortedSet, i) = ((start, stop) = get_range(rs.ranges, i); stop - start + 1)
@inline get_range_indices(rs::UnitRangesSortedVector, i) = (r = rs.ranges[i]; return (r.start, r.stop))
@inline get_range_indices(rs::UnitRangesSortedSet, i::DataStructures.Tokens.IntSemiToken) =
    ((start, stop) = deref((rs.ranges, i)); return (start, stop))
@inline get_range(rs::UnitRangesSortedVector, i) = rs.ranges[i]
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
@inline get_range_start(rs::UnitRangesSortedVector, i) = rs.ranges[i].start
@inline get_range_start(rs::UnitRangesSortedSet, i) = deref_key((rs.ranges, i))
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

@inline Base.firstindex(rs::UnitRangesSortedVector) = firstindex(rs.ranges)
@inline Base.firstindex(rs::UnitRangesSortedSet) = startof(rs.ranges)
@inline Base.lastindex(rs::AbstractUnitRangesSortedSet) = lastindex(rs.ranges)
@inline Base.first(rs::UnitRangesSortedVector) = rs.ranges[1]
@inline Base.first(rs::UnitRangesSortedSet) = ((start,stop) = deref((rs.ranges, startof(rs.ranges))); (start:stop))
@inline Base.last(rs::UnitRangesSortedVector) = rs.ranges[end]
@inline Base.last(rs::UnitRangesSortedSet) = ((start,stop) = deref((rs.ranges, lastindex(rs.ranges))); (start:stop))

@inline beforestartindex(rs::UnitRangesSortedVector) = firstindex(rs.ranges) - 1
@inline beforestartindex(rs::UnitRangesSortedSet) = beforestartsemitoken(rs.ranges)
@inline pastendindex(rs::UnitRangesSortedVector) = lastindex(rs.ranges) + 1
@inline pastendindex(rs::UnitRangesSortedSet) = pastendsemitoken(rs.ranges)

@inline DataStructures.advance(rs::UnitRangesSortedVector, state) = state + 1
@inline DataStructures.advance(rs::UnitRangesSortedSet, state) = advance((rs.ranges, state))

@inline urangevectorisless(x::UnitRange, y::UnitRange) = x.start < y.start
@inline urangevectorisless(x::UnitRange, y) = x.start < y
@inline urangevectorisless(x, y::UnitRange) = x < y.start

@inline searchsortedlastrange(rs::UnitRangesSortedVector, i) = searchsortedlast(rs.ranges, i; lt=urangevectorisless)
@inline searchsortedlastrange(rs::UnitRangesSortedSet, i) = searchsortedlast(rs.ranges, i)

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


Base.@propagate_inbounds function check_exist_and_update(rs, i)

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
        push!(rs.ranges, (i:i))
        rs.lastusedrangeindex = 1
        return rs
    end

    if st == beforestartindex(rs)  # the index `i` is before the first range
        if rs.ranges[1].start - i > 1  # there is will be gap in indices after inserting
            pushfirst!(rs.ranges, (i:i))
        else  # prepend to first range
            rs.ranges[1] = (i : rs.ranges[1].stop)
        end
        rs.lastusedrangeindex = 1
        return rs
    end

    r = get_range(rs, st)

    if i >= rs.ranges[end].start  # the index `i` is after the last range start
        if i > r.stop + 1  # there is will be the gap in indices after inserting
            push!(rs.ranges, (i:i))
        else  # just append to last range
            rs.ranges[st] = (r.start : r.stop+1)
        end
        rs.lastusedrangeindex = length(rs.ranges)
        return rs
    end

    # the index `i` is somewhere between indices
    stnext = advance(rs, st)
    rnext = rs.ranges[stnext]

    if rnext.start - r.stop == 2  # join ranges
        rs.ranges[st] = (r.start:rnext.stop)
        deleteat!(rs.ranges, stnext)
        rs.lastusedrangeindex = st
    elseif i - r.stop == 1  # append to left range
        rs.ranges[st] = (r.start : r.stop+1)
        rs.lastusedrangeindex = st
    elseif rnext.start - i == 1  # prepend to right range
        rs.ranges[stnext] = (rnext.start-1 : rnext.stop)
        rs.lastusedrangeindex = stnext
    else  # insert single element chunk
        insert!(rs.ranges, stnext, (i:i))
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



@inline function Base.delete!(rs::UnitRangesSortedVector{Ti}, idx) where {Ti}
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
        deleteat!(rs.ranges, st)
    elseif i == r.stop  # last index in range
        rs.ranges[st] = (r.start : r.stop-1)
    elseif i == r.start  # first index in range
        rs.ranges[st] = (r.start+1 : r.stop)
    else
        insert!(rs.ranges, advance(rs, st), (i+1 : r.stop))
        rs.ranges[st] = (r.start : i-1)
    end

    rs.lastusedrangeindex = beforestartindex(rs)

    return rs
end

@inline function Base.delete!(rs::UnitRangesSortedSet{Ti}, idx) where {Ti}
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


function Base.empty!(rs::AbstractUnitRangesSortedSet)
    empty!(rs.ranges)
    rs.lastusedrangeindex = beforestartindex(rs)
    rs
end


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
