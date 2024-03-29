
using BenchmarkTools
using DataStructures
using UnitRangesSortedSets
using UnitRangesSortedSets: advance, beforefirstindex, inferrangetype, pastlastindex, regress, to_urange
using Test

using Aqua
Aqua.test_all(UnitRangesSortedSets)

# TODO:
# Custom Number Type, see Unitful: "runtests.jl"
#
# Try Test.GenericSet
#
# [rs[i] for i in searchsortedrange(rs, 10:40)]
#
# function VecUnitRangesSortedSet(rs::AbstractUnitRangesSortedSet{K,TU}) where {K,TU}
# function UnitRangesSortedSet(rs::AbstractUnitRangesSortedSet{K,TU}) where {K,TU}
# @inline function URSSIndexURange(rs::P, l::Tix, r::Tix) where {P<:VecUnitRangesSortedSet, Tix}
# URSSIndexURange for VecUnitRangesSortedSet?
# @inline Base.firstindex(ur::URSSIndexURange) = ur.start
# @inline Base.lastindex(ur::URSSIndexURange) = ur.stop
# @inline DataStructures.advance(ur::URSSIndexURange), beforefirstindex
# @inline function Base.iterate(ur::URSSIndexURange, state = (first(ur), 0))
# @inline function Base.getindex(ur::URSSIndexURange, i) ???
# @inline function indexcompare(rs::VecUnitRangesSortedSet, i, j)
# @inline function indexcompare(rs::AbstractSubUnitRangesSortedSet{K,TU,P}, i, j) where {K,TU,P<:VecUnitRangesSortedSet}
# @inline getindex_tuple(rs::Sub01XUnitRangesSortedSet, i)
# @inline searchsortedrangelast(rs::Sub01XUnitRangesSortedSet, k)
# @inline searchsortedrangefirst(rs::Sub01XUnitRangesSortedSet, k)
# @inline searchsortedrange(rs::AbstractSubUnitRangesSortedSet, kk::AbstractRange)
#
# function Base.show(io::IO, x::T) where {T<:AbstractUnitRangesSortedSubSet}
# function Base.show(io::IO, ::MIME"text/plain", x::T) where {T<:AbstractUnitRangesSortedSubSet}
#



# https://github.com/JuliaLang/julia/issues/39952
basetype(::Type{T}) where T = Base.typename(T).wrapper


const list_of_Ti_to_test = (Int, UInt16, Float64, Char)
#const list_of_Ti_to_test = (Int, UInt64, UInt16, Float64, Char)

const list_of_containers_types_to_test = (VecUnitRangesSortedSet, UnitRangesSortedSet)


function check_ranges_isequal(r, rs...)
    for rr in rs
        step(r)  == step(rr)  || return false
        first(r) == first(rr) || return false
        last(r)  == last(rr)  || return false
    end
    return true
end

function check_isequal(rs1, rs2)
    length(rs1) == length(rs2) || return false
    for (r1, r2) in zip(rs1, rs2)
        check_ranges_isequal(r1, r2) || return false
    end
    return true
end

function test_isequal(rs1, rs2)
    @test length(rs1) == length(rs2)
    for (r1, r2) in zip(rs1, rs2)
        @test step(r1)  == step(r2)
        @test first(r1) == first(r2)
        @test last(r1)  == last(r2)
    end
end

function test_range_create(rs::T, vu, Tv::Type) where {T<:AbstractUnitRangesSortedContainer}
    TRange = inferrangetype(Tv)
    @test eltype(rs) === TRange
    @test typeof(rs) === basetype(T){Tv,TRange}
    test_isequal(rs, vu)
end


convertinfer(K::Type, rs::NTuple{N,T}) where{N,T<:AbstractRange} =
    convert(NTuple{length(rs),inferrangetype(K)}, rs)

convertinfer(K::Type, rs::NTuple{N,T}) where{N,T} =
    convert(NTuple{length(rs),K}, rs)

convertinfer(K::Type, rs::AbstractVector{T}) where{T<:AbstractRange} =
    convert(Vector{inferrangetype(K)}, rs)

convertinfer(K::Type, rs::AbstractVector{T}) where{T} =
    convert(Vector{K}, rs)

convertinfer(K::Type, r::AbstractRange) =
    to_urange(inferrangetype(K), r)



function test_iterators(rs, vu)
    @test length(rs) == length(vu)
    @test regress(rs, firstindex(rs)) == beforefirstindex(rs)
    @test advance(rs, lastindex(rs)) == pastlastindex(rs)
    @test check_ranges_isequal(rs[firstindex(rs)], rs[begin], first(rs), first(vu))
    @test check_ranges_isequal(rs[lastindex(rs)], rs[end], last(rs), last(vu))
    @test check_isequal(rs, vu)
    @test check_isequal(Iterators.reverse(rs), Iterators.reverse(vu))
    @test check_isequal([r for r in rs], [r for r in vu])
    @test check_isequal([r for r in Iterators.reverse(rs)], [r for r in Iterators.reverse(vu)])
    @test sum(r->length(r), rs) ==
          sum(r->length(r), Iterators.reverse(rs)) ==
          sum(i->length(rs[i]), eachindex(rs)) ==
          sum(i->length(rs[i]), eachindex(Iterators.reverse(rs))) ==
          sum(i->length(rs[i]), Iterators.reverse(eachindex(rs))) ==
          sum(r->length(r), vu)
end

function test_empty_iterators(rs, vu)
    @test length(rs) == length(vu)
    @test firstindex(rs) == pastlastindex(rs)
    @test lastindex(rs) == beforefirstindex(rs)
    @test_throws BoundsError rs[begin]
    @test_throws BoundsError rs[end]
    @test_throws BoundsError first(rs)
    @test_throws BoundsError last(rs)
    @test check_isequal(rs, vu)
    @test check_isequal(Iterators.reverse(rs), Iterators.reverse(vu))
    @test check_isequal([r for r in rs], [r for r in vu])
    @test check_isequal([r for r in Iterators.reverse(rs)], [r for r in Iterators.reverse(vu)])
    @test sum(r->length(r), rs, init=0) ==
          sum(r->length(r), Iterators.reverse(rs), init=0) ==
          sum(i->length(rs[i]), eachindex(rs), init=0) ==
          sum(i->length(rs[i]), eachindex(Iterators.reverse(rs)), init=0) ==
          sum(i->length(rs[i]), Iterators.reverse(eachindex(rs)), init=0) ==
          sum(r->length(r), vu, init=0)
end



@testset "Creating" begin
    for K in list_of_Ti_to_test
        for TypeURSS in list_of_containers_types_to_test
            @eval begin

                vu = Vector{inferrangetype($K)}(undef, 0)
                test_range_create($TypeURSS{$K}(), convertinfer($K, vu), $K)

                v = [1, 2]
                vu = [1:2]
                test_range_create($TypeURSS{$K}(v), convertinfer($K, vu), $K)

                v = [1, 3, 4]
                vu = [1:1, 3:4]
                test_range_create($TypeURSS{$K}(v), convertinfer($K, vu), $K)

                v = [3:4, 1:1]
                vu = [1:1, 3:4]
                test_range_create($TypeURSS{$K}(v), convertinfer($K, vu), $K)

                v = [3:4, 2:5, 1:1]
                vu = [1:5]
                test_range_create($TypeURSS{$K}(v), convertinfer($K, vu), $K)

                v = [$K(1), $K(3), $K(4)]
                vu = [1:1, 3:4]
                test_range_create($TypeURSS(v), convertinfer($K, vu), $K)

                v = [$K(3):$K(4), $K(1):$K(1)]
                vu = [1:1, 3:4]
                test_range_create($TypeURSS(v), convertinfer($K, vu), $K)

                v = [UInt64(3):UInt64(4), UInt64(1):UInt64(1)]
                vu = [1:1, 3:4]
                test_range_create($TypeURSS{$K}(v), convertinfer($K, vu), $K)

                v = [UInt64(3), UInt64(4), UInt64(1)]
                vu = [1:1, 3:4]
                test_range_create($TypeURSS(v), convertinfer(UInt64, vu), UInt64)

                v = [UInt64(3):UInt64(4), UInt64(1):UInt64(1)]
                vu = [1:1, 3:4]
                test_range_create($TypeURSS(v), convertinfer(UInt64, vu), UInt64)
            end
        end
    end
end


@testset "Searching" begin
    for K in list_of_Ti_to_test
        for TypeURSS in list_of_containers_types_to_test
            @eval begin
                rs = $TypeURSS{$K}((1:2, 4:4, 6:6))

                ir = searchsortedrangelast(rs, 0)
                @test ir == beforefirstindex(rs)

                ir = searchsortedrangelast(rs, 1)
                @test ir != beforefirstindex(rs) &&
                      getindex(rs, ir) == convertinfer($K, 1:2) &&
                      ir == firstindex(rs)

                ir = searchsortedrangelast(rs, 2)
                @test ir != beforefirstindex(rs) &&
                      getindex(rs, ir) == convertinfer($K, 1:2) &&
                      ir == firstindex(rs)

                ir = searchsortedrangelast(rs, 3)
                @test ir != beforefirstindex(rs) &&
                      getindex(rs, ir) == convertinfer($K, 1:2) &&
                      ir == firstindex(rs)

                ir = searchsortedrangelast(rs, 4)
                @test ir != beforefirstindex(rs) &&
                      getindex(rs, ir) == convertinfer($K, 4:4) &&
                      ir != firstindex(rs) && ir != lastindex(rs)

                ir = searchsortedrangelast(rs, 5)
                @test ir != beforefirstindex(rs) &&
                      getindex(rs, ir) == convertinfer($K, 4:4) &&
                      ir != firstindex(rs) && ir != lastindex(rs)

                ir = searchsortedrangelast(rs, 6)
                @test ir != beforefirstindex(rs) &&
                      getindex(rs, ir) == convertinfer($K, 6:6) &&
                      ir == lastindex(rs)

                ir = searchsortedrangelast(rs, 7)
                @test ir != beforefirstindex(rs) &&
                      getindex(rs, ir) == convertinfer($K, 6:6) &&
                      ir == lastindex(rs)


                ir = searchsortedrangefirst(rs, 0)
                @test ir != pastlastindex(rs) &&
                      getindex(rs, ir) == convertinfer($K, 1:2) &&
                      ir == firstindex(rs)

                ir = searchsortedrangefirst(rs, 1)
                @test ir != pastlastindex(rs) &&
                      getindex(rs, ir) == convertinfer($K, 1:2) &&
                      ir == firstindex(rs)

                ir = searchsortedrangefirst(rs, 2)
                @test ir != pastlastindex(rs) &&
                      getindex(rs, ir) == convertinfer($K, 1:2) &&
                      ir == firstindex(rs)

                ir = searchsortedrangefirst(rs, 3)
                @test ir != pastlastindex(rs) &&
                      getindex(rs, ir) == convertinfer($K, 4:4) &&
                      ir != firstindex(rs) && ir != lastindex(rs)

                ir = searchsortedrangefirst(rs, 4)
                @test ir != pastlastindex(rs) &&
                      getindex(rs, ir) == convertinfer($K, 4:4) &&
                      ir != firstindex(rs) && ir != lastindex(rs)

                ir = searchsortedrangefirst(rs, 5)
                @test ir != pastlastindex(rs) &&
                      getindex(rs, ir) == convertinfer($K, 6:6) &&
                      ir == lastindex(rs)

                ir = searchsortedrangefirst(rs, 6)
                @test ir != pastlastindex(rs) &&
                      getindex(rs, ir) == convertinfer($K, 6:6) &&
                      ir == lastindex(rs)

                ir = searchsortedrangefirst(rs, 7)
                @test ir == pastlastindex(rs)



                rs = $TypeURSS{$K}((1:2, 4:4))

                iir = searchsortedrange(rs, 0)
                @test length(iir) == 0 && getindex(rs, first(iir)) == convertinfer($K, 1:2) &&
                                         first(iir) == firstindex(rs) &&
                                         last(iir) == beforefirstindex(rs)

                iir = searchsortedrange(rs, 1)
                @test length(iir) == 1 && getindex(rs, first(iir)) == convertinfer($K, 1:2) &&
                                         first(iir) == firstindex(rs)

                iir = searchsortedrange(rs, 2)
                @test length(iir) == 1 && getindex(rs, first(iir)) == convertinfer($K, 1:2)

                iir = searchsortedrange(rs, 3)
                @test length(iir) == 0 && getindex(rs, first(iir)) == convertinfer($K, 4:4) &&
                                         getindex(rs, last(iir)) == convertinfer($K, 1:2)

                iir = searchsortedrange(rs, 4)
                @test length(iir) == 1 && getindex(rs, first(iir)) == convertinfer($K, 4:4)

                iir = searchsortedrange(rs, 5)
                @test length(iir) == 0 && getindex(rs, last(iir)) == convertinfer($K, 4:4) &&
                                         first(iir) == pastlastindex(rs) &&
                                         last(iir) == lastindex(rs)

                rs = $TypeURSS{$K}((1:2, 4:4))
                @test getrange(rs, 0) === nothing
                @test getrange(rs, 1) == convertinfer($K, 1:1)
                @test getrange(rs, 1:1) == convertinfer($K, 1:1)
                @test getrange(rs, 1:2) == convertinfer($K, 1:2)
                @test getrange(rs, 1:3) === nothing
            end
        end
    end
end


@testset "`in`" begin
    for K in list_of_Ti_to_test
        for TypeURSS in list_of_containers_types_to_test
            @eval begin

                v = [1, 2, 3, 5, 6, 8]
                rs = $TypeURSS{$K}(v)
                for i = 0:9
                    @test !xor(in(i, rs), in(i, v))
                end

                v = [2:3, 5:6]
                rs = $TypeURSS{$K}(v)
                @test !in(1:2, rs)
                @test in(2:2, rs)
                @test in(2:3, rs)
                @test in(5:5, rs)
                @test !in(2:4, rs)
                @test !in(4:4, rs)
                @test !in(2:6, rs)
                @test !in(1:7, rs)
            end
        end
    end
end


@testset "`push!`" begin
    for K in list_of_Ti_to_test
        for TypeURSS in list_of_containers_types_to_test
            @eval begin

                rs = $TypeURSS{$K}()
                @test check_isequal(push!(rs, 3),  convertinfer($K, (3:3,)))
                @test check_isequal(push!(rs, 4),  convertinfer($K, (3:4,)))
                @test check_isequal(push!(rs, 18), convertinfer($K, (3:4, 18:18)))
                @test check_isequal(push!(rs, 2),  convertinfer($K, (2:4, 18:18)))
                @test check_isequal(push!(rs, 0),  convertinfer($K, (0:0, 2:4, 18:18)))
                @test check_isequal(push!(rs, 20), convertinfer($K, (0:0, 2:4, 18:18, 20:20)))
                @test check_isequal(push!(rs, 19), convertinfer($K, (0:0, 2:4, 18:20)))
                @test check_isequal(push!(rs, 1),  convertinfer($K, (0:4, 18:20)))
                @test check_isequal(push!(rs, 9),  convertinfer($K, (0:4, 9:9, 18:20)))
                @test check_isequal(push!(rs, 9),  convertinfer($K, (0:4, 9:9, 18:20)))
                @test check_isequal(push!(rs, 10), convertinfer($K, (0:4, 9:10, 18:20)))
                @test check_isequal(push!(rs, 8),  convertinfer($K, (0:4, 8:10, 18:20)))

                rs = $TypeURSS{$K}()
                @test check_isequal(push!(rs, 13:14), convertinfer($K, (13:14,)))
                @test check_isequal(push!(rs, 23:24), convertinfer($K, (13:14, 23:24)))
                @test check_isequal(push!(rs, 12:12), convertinfer($K, (12:14, 23:24)))
                @test check_isequal(push!(rs, 11:13), convertinfer($K, (11:14, 23:24)))
                @test check_isequal(push!(rs, 8:9),   convertinfer($K, (8:9, 11:14, 23:24)))
                @test check_isequal(push!(copy(rs), [9:12]),  convertinfer($K, (8:14, 23:24)))
                @test check_isequal(push!(copy(rs), (9:12,)),  convertinfer($K, (8:14, 23:24)))
                @test check_isequal(push!(copy(rs), [9,10,11,12]),  convertinfer($K, (8:14, 23:24)))
                @test check_isequal(push!(rs, 9:12),  convertinfer($K, (8:14, 23:24)))
                @test check_isequal(push!(rs, 16:17),  convertinfer($K, (8:14, 16:17, 23:24)))
                @test check_isequal(push!(rs, 25:26),  convertinfer($K, (8:14, 16:17, 23:26)))
                @test check_isequal(push!(rs, 24:26),  convertinfer($K, (8:14, 16:17, 23:26)))
                @test check_isequal(push!(rs, 28:29),  convertinfer($K, (8:14, 16:17, 23:26, 28:29)))
                @test check_isequal(push!(rs, 28:29),  convertinfer($K, (8:14, 16:17, 23:26, 28:29)))
                @test check_isequal(push!(rs, 16:15),  convertinfer($K, (8:14, 16:17, 23:26, 28:29)))
                @test check_isequal(push!(rs, 18:19),  convertinfer($K, (8:14, 16:19, 23:26, 28:29)))
            end
        end
    end
end


@testset "`delete!`" begin
    for K in list_of_Ti_to_test
        for TypeURSS in list_of_containers_types_to_test
            @eval begin

                rs = $TypeURSS{$K}((0:0, 2:4, 8:8, 10:12))
                @test check_isequal(delete!(rs, 5),  convertinfer($K, (0:0, 2:4, 8:8, 10:12)))
                @test check_isequal(delete!(rs, 0),  convertinfer($K, (2:4, 8:8, 10:12)))
                @test check_isequal(delete!(rs, 0),  convertinfer($K, (2:4, 8:8, 10:12)))
                @test check_isequal(delete!(rs, 8),  convertinfer($K, (2:4, 10:12)))
                @test check_isequal(delete!(rs, 12), convertinfer($K, (2:4, 10:11)))
                @test check_isequal(delete!(rs, 10), convertinfer($K, (2:4, 11:11)))
                @test check_isequal(delete!(rs, 11), convertinfer($K, (2:4,)))
                @test check_isequal(delete!(rs, 3),  convertinfer($K, (2:2, 4:4)))
                @test delete!(rs, 3) === rs

                rs = $TypeURSS{$K}((0:0, 2:6, 8:8, 10:13, 15:20))
                @test check_isequal(delete!(rs, 17:17),  convertinfer($K, (0:0, 2:6, 8:8, 10:13, 15:16, 18:20)))
                @test check_isequal(delete!(rs, 16:18),  convertinfer($K, (0:0, 2:6, 8:8, 10:13, 15:15, 19:20)))
                @test check_isequal(delete!(rs, 15:20),  convertinfer($K, (0:0, 2:6, 8:8, 10:13)))
                @test check_isequal(delete!(rs, 7:7),  convertinfer($K, (0:0, 2:6, 8:8, 10:13)))
                @test check_isequal(delete!(rs, 6:7),  convertinfer($K, (0:0, 2:5, 8:8, 10:13)))
                @test check_isequal(delete!(rs, 0:1),  convertinfer($K, (2:5, 8:8, 10:13)))
                @test check_isequal(delete!(rs, 9:10),  convertinfer($K, (2:5, 8:8, 11:13)))
                @test check_isequal(delete!(rs, 13:14),  convertinfer($K, (2:5, 8:8, 11:12)))
                @test check_isequal(delete!(copy(rs), [1:3]),  convertinfer($K, (4:5, 8:8, 11:12)))
                @test check_isequal(delete!(copy(rs), (1:3,)),  convertinfer($K, (4:5, 8:8, 11:12)))
                @test check_isequal(delete!(copy(rs), [1,2,3]),  convertinfer($K, (4:5, 8:8, 11:12)))
                @test check_isequal(delete!(copy(rs), (1,2,3)),  convertinfer($K, (4:5, 8:8, 11:12)))
                @test check_isequal(delete!(rs, 1:3),  convertinfer($K, (4:5, 8:8, 11:12)))

            end
        end
    end
end


@testset "`subset`" begin
    for K in list_of_Ti_to_test
        for TypeURSS in list_of_containers_types_to_test
            @eval begin

                rs = $TypeURSS{$K}((0:0, 2:6, 8:16, 20:33, 35:47, 49:50))

                it = searchsortedrangelast(rs, 1)

                ss = subset(rs, :)
                @test typeof(ss) == UnitRangesSortedSubSet{$K,inferrangetype($K),typeof(rs),typeof(it)}
                vu = (0:0, 2:6, 8:16, 20:33, 35:47, 49:50)
                test_iterators(ss, convertinfer($K, vu))

                ss = subset(rs, 0:50)
                vu = (0:0, 2:6, 8:16, 20:33, 35:47, 49:50)
                test_iterators(ss, convertinfer($K, vu))

                ss = subset(rs, 0:0)
                @test typeof(ss) == UnitRangesSortedSubSet1{$K,inferrangetype($K),typeof(rs),typeof(it)}
                vu = (0:0,)
                test_iterators(ss, convertinfer($K, vu))

                ss = subset(rs, 0:1)
                @test typeof(ss) == UnitRangesSortedSubSet1{$K,inferrangetype($K),typeof(rs),typeof(it)}
                vu = (0:0,)
                test_iterators(ss, convertinfer($K, vu))

                ss = subset(rs, 10:40)
                @test typeof(ss) == UnitRangesSortedSubSet{$K,inferrangetype($K),typeof(rs),typeof(it)}
                vu = (10:16, 20:33, 35:40)
                test_iterators(ss, convertinfer($K, vu))

                ss = subset(rs, 7:40)
                vu = (8:16, 20:33, 35:40)
                test_iterators(ss, convertinfer($K, vu))

                ss = subset(rs, 11:48)
                vu = (11:16, 20:33, 35:47)
                test_iterators(ss, convertinfer($K, vu))

                ss = subset(rs, 7:48)
                vu = (8:16, 20:33, 35:47)
                test_iterators(ss, convertinfer($K, vu))

                ss = subset(rs, 21:32)
                vu = (21:32,)
                test_iterators(ss, convertinfer($K, vu))

                ss = subset(rs, 20:33)
                vu = (20:33,)
                test_iterators(ss, convertinfer($K, vu))

                ss = subset(rs, 19:34)
                vu = (20:33,)
                test_iterators(ss, convertinfer($K, vu))

                ss = subset(rs, 22:34)
                vu = (22:33,)
                test_iterators(ss, convertinfer($K, vu))

                ss = subset(rs, 19:30)
                vu = (20:30,)
                test_iterators(ss, convertinfer($K, vu))

                ss = subset(rs, 18:19)
                @test typeof(ss) == UnitRangesSortedSubSet0{$K,inferrangetype($K),typeof(rs),typeof(it)}
                vu = Vector{eltype(ss)}([])
                test_empty_iterators(ss, convertinfer($K, vu))

                ss = subset(rs, 21:20)
                @test typeof(ss) == UnitRangesSortedSubSet0{$K,inferrangetype($K),typeof(rs),typeof(it)}
                vu = Vector{eltype(ss)}([])
                test_empty_iterators(ss, convertinfer($K, vu))

                ss = subset(rs, 51:51)
                @test typeof(ss) == UnitRangesSortedSubSet0{$K,inferrangetype($K),typeof(rs),typeof(it)}
                vu = Vector{eltype(ss)}([])
                test_empty_iterators(ss, convertinfer($K, vu))

                delete!(rs, 0:0)
                ss = subset(rs, 0:1)
                @test typeof(ss) == UnitRangesSortedSubSet0{$K,inferrangetype($K),typeof(rs),typeof(it)}
                vu = Vector{eltype(ss)}([])
                test_empty_iterators(ss, convertinfer($K, vu))

                # subset of subset
                rs = $TypeURSS{$K}((0:0, 2:6, 8:16, 20:33, 35:47, 49:50))

                it = searchsortedrangelast(rs, 1)

                ss = subset(rs, :)

                sss = subset(ss, :)
                @test typeof(sss) == UnitRangesSortedSubSet{$K,inferrangetype($K),typeof(rs),typeof(it)}
                vu = (0:0, 2:6, 8:16, 20:33, 35:47, 49:50)
                test_iterators(sss, convertinfer($K, vu))

                sss = subset(ss, 0:0)
                @test typeof(sss) == UnitRangesSortedSubSet1{$K,inferrangetype($K),typeof(rs),typeof(it)}
                vu = (0:0,)
                test_iterators(sss, convertinfer($K, vu))

                ss = subset(rs, 10:40) # (10:16, 20:33, 35:40)

                sss = subset(ss, 15:35)
                @test typeof(sss) == UnitRangesSortedSubSet{$K,inferrangetype($K),typeof(rs),typeof(it)}
                vu = (15:16, 20:33, 35:35)
                test_iterators(sss, convertinfer($K, vu))

                sss = subset(ss, 18:19)
                @test typeof(sss) == UnitRangesSortedSubSet0{$K,inferrangetype($K),typeof(rs),typeof(it)}
                vu = Vector{eltype(ss)}([])
                test_empty_iterators(sss, convertinfer($K, vu))

            end
        end
    end
end


@testset "`in` subsets" begin
    for K in list_of_Ti_to_test
        for TypeURSS in list_of_containers_types_to_test
            @eval begin

                v = [0, 1, 2, 3, 5, 6, 8, 9, 11, 12]
                rs = $TypeURSS{$K}(v)
                vu = [1, 2, 3, 5, 6, 8]
                ss = subset(rs, 1:8)
                for i = 0:13
                    @test !xor(in(i, ss), in(i, vu))
                end

                v = [0:0, 2:3, 5:6, 8:8]
                rs = $TypeURSS{$K}(v)
                ss = subset(rs, 2:6)
                @test !in(1:2, ss)
                @test in(2:2, ss)
                @test in(2:3, ss)
                @test in(5:5, ss)
                @test !in(2:4, ss)
                @test !in(4:4, ss)
                @test !in(2:6, ss)
                @test !in(1:7, ss)

                v = [0:0, 2:3, 5:6, 8:8]
                rs = $TypeURSS{$K}(v)
                ss = subset(rs, 2:3)
                @test !in(1:2, ss)
                @test !in(1, ss)
                @test in(2, ss)
                @test in(2:2, ss)
                @test in(2:3, ss)
                @test in(3:3, ss)
                @test !in(2:4, ss)
                @test !in(3:4, ss)
                @test !in(4:4, ss)
                @test !in(5:5, ss)

                v = [0:0, 2:3, 5:6, 8:8]
                rs = $TypeURSS{$K}(v)
                ss = subset(rs, 3:2)
                @test !in(1, ss)
                @test !in(2, ss)
                @test !in(1:2, ss)
                @test !in(2:2, ss)
                @test !in(2:3, ss)
                @test !in(3:3, ss)
                @test !in(2:4, ss)
                @test !in(3:4, ss)
                @test !in(4:4, ss)
                @test !in(5:5, ss)
            end
        end
    end
end


@testset "Common functions" begin
    for K in list_of_Ti_to_test
        for TypeURSS in list_of_containers_types_to_test
            @eval begin

                rs = $TypeURSS{$K}()
                vu = Vector{eltype(rs)}([])
                test_empty_iterators(rs, convertinfer($K, vu))

                rs = $TypeURSS{$K}((1:6,))
                vu = (1:6,)
                test_iterators(rs, convertinfer($K, vu))

                rs = $TypeURSS{$K}((1:6, 8:16, 20:33))
                vu = (1:6, 8:16, 20:33)
                test_iterators(rs, convertinfer($K, vu))


                rs = $TypeURSS{$K}((0:0, 2:4))
                rs2 = copy(rs)
                @test rs2 == rs && rs2 !== rs
                rs2 = empty(rs)
                @test rs2 == $TypeURSS{$K}()
                empty!(rs2)
                @test rs2 != rs && rs2 !== rs
                @test length(rs2) == 0

                rs = $TypeURSS{$K}((0:0, 2:4))
                @test filter(s->length(s)==3, rs) == $TypeURSS{$K}((2:4))
                @test (filter!(s->length(s)==1, rs); rs == $TypeURSS{$K}((0:0)))

                rs = $TypeURSS{$K}((0:0, 2:4))
                @test findfirst(s->issubset(s, $K(0):$K(4)), rs) == searchsortedrangelast(rs, $K(0))
                @test findfirst(s->issubset(s, $K(1):$K(4)), rs) == searchsortedrangelast(rs, $K(2))
                @test findall(s->length(s)==1, rs) == convertinfer($K, [0:0])

                v = [0:0, 2:4]
                rs = $TypeURSS{$K}(v)
                @test collect(rs) == convertinfer($K, v)
                @test collect(UInt64, rs) == convertinfer(UInt64, v)

                @test convert(typeof(rs), rs) === rs

                @test convert(SortedSet{inferrangetype($K)}, rs) == SortedSet(to_urange.(inferrangetype($K), v))
                @test convert(SortedSet{UInt64}, rs)             == SortedSet{UInt64}(collect(Iterators.flatten(v)))
                @test convert(SortedSet{$K}, rs)                 == SortedSet{$K}(collect(Iterators.flatten(v)))
                @test convert(SortedSet, rs)                     == SortedSet{$K}(collect(Iterators.flatten(v)))
                @test convert(Vector{inferrangetype($K)}, rs)    == Vector(to_urange.(inferrangetype($K), v))
                @test convert(Vector{UInt64}, rs)                == Vector{UInt64}(collect(Iterators.flatten(v)))
                @test convert(Vector{$K}, rs)                    == Vector{$K}(collect(Iterators.flatten(v)))
                @test convert(Vector, rs)                        == Vector{$K}(collect(Iterators.flatten(v)))
            end
        end
    end
end


@testset "`Set` API" begin
    for K in list_of_Ti_to_test
        for TypeURSS in list_of_containers_types_to_test
            @eval begin

                @test union($TypeURSS{$K}((0:0, 2:4)), $TypeURSS{$K}((2:3, 5:6))) == $TypeURSS{$K}((0:0, 2:6))
                @test union($TypeURSS{$K}((0:0, 2:4)), (2:3, 5:6)) == $TypeURSS{$K}((0:0, 2:6))
                @test union($TypeURSS{$K}((0:0, 2:4)), (2, 3, 5, 6)) == $TypeURSS{$K}((0:0, 2:6))
                @test union($TypeURSS{$K}((0:0, 2:4)), [2:3], [5:6]) == $TypeURSS{$K}((0:0, 2:6))
                @test union($TypeURSS{$K}((0:0, 2:4)), [2:3, 5:6]) == $TypeURSS{$K}((0:0, 2:6))
                @test union($TypeURSS{$K}((0:0, 2:4)), [2, 3, 5, 6]) == $TypeURSS{$K}((0:0, 2:6))
                @test union($TypeURSS{$K}((0:0, 2:4)), 1:1) == $TypeURSS{$K}((0:4))
                @test union($TypeURSS{$K}((0:0, 2:4)), 1) == $TypeURSS{$K}((0:4))

                @test intersect($TypeURSS{$K}((0:0, 2:4)), $TypeURSS{$K}((2:3, 5:6))) == $TypeURSS{$K}((2:3))
                @test intersect($TypeURSS{$K}((0:0, 2:4)), (2:3, 5:6)) == $TypeURSS{$K}((2:3))
                @test intersect($TypeURSS{$K}((0:0, 2:4)), 2:3) == $TypeURSS{$K}((2:3))
                @test intersect($TypeURSS{$K}((0:0, 2:4)), 2) == $TypeURSS{$K}((2:2))
                rs = $TypeURSS{$K}((0:0, 2:4))
                @test (intersect!(rs, $TypeURSS{$K}((2:3, 5:6, 8:8))); rs == $TypeURSS{$K}((2:3)))
                rs = [0:0, 2:4, 8:10, 12:12]
                @test (intersect(rs, $TypeURSS{$K}((2:3, 5:6, 8:8))) == [2:3, 8:8])
                @test (intersect!(rs, $TypeURSS{$K}((2:3, 5:6, 8:8))); rs == [2:3, 8:8])
                rs = Set((0:0, 2:4, 8:10, 12:12))
                @test (intersect(rs, $TypeURSS{$K}((2:3, 5:6, 8:8))) == Set((2:3, 8:8)))
                @test (intersect!(rs, $TypeURSS{$K}((2:3, 5:6, 8:8))); rs == Set((2:3, 8:8)))
                rs = SortedSet((0:0, 2:4, 8:10, 12:12))
                @test (intersect(rs, $TypeURSS{$K}((2:3, 5:6, 8:8))) == SortedSet((2:3, 8:8)))
                @test (intersect!(rs, $TypeURSS{$K}((2:3, 5:6, 8:8))); rs == SortedSet((2:3, 8:8)))
                rs = [0, 2,3,4, 8,9,10, 12]
                @test (intersect(rs, $TypeURSS{$K}((2:3, 5:6, 8:8))) == [2,3, 8])
                @test (intersect!(rs, $TypeURSS{$K}((2:3, 5:6, 8:8))); rs == [2,3, 8])
                rs = Set((0, 2,3,4, 8,9,10, 12))
                @test (intersect(rs, $TypeURSS{$K}((2:3, 5:6, 8:8))) == Set((2,3, 8)))
                @test (intersect!(rs, $TypeURSS{$K}((2:3, 5:6, 8:8))); rs == Set((2,3, 8)))
                rs = SortedSet((0, 2,3,4, 8,9,10, 12))
                @test (intersect(rs, $TypeURSS{$K}((2:3, 5:6, 8:8))) == SortedSet((2,3, 8)))
                @test (intersect!(rs, $TypeURSS{$K}((2:3, 5:6, 8:8))); rs == SortedSet((2,3, 8)))

                rs = $TypeURSS{$K}((2:3, 5:6, 8:8))
                @test (intersect!(rs, convertinfer($K, [0:0, 2:4, 8:10, 12:12])); rs == $TypeURSS{$K}((2:3, 8:8)))
                rs = $TypeURSS{$K}((2:3, 5:6, 8:8))
                @test (intersect!(rs, convertinfer($K, [0, 2,3,4, 8,9,10, 12])); rs == $TypeURSS{$K}((2:3, 8:8)))
                rs = $TypeURSS{$K}((2:3, 5:6, 8:8))
                @test (intersect!(rs, Set(convertinfer($K, [0:0, 2:4, 8:10, 12:12]))); rs == $TypeURSS{$K}((2:3, 8:8)))
                rs = $TypeURSS{$K}((2:3, 5:6, 8:8))
                @test (intersect!(rs, Set(convertinfer($K, [0, 2,3,4, 8,9,10, 12]))); rs == $TypeURSS{$K}((2:3, 8:8)))


                rs = $TypeURSS{$K}((2:3, 5:6))
                @test (setdiff(rs, Set(2:2)) == $TypeURSS{$K}((3:3, 5:6)))
                @test (setdiff(rs, OrderedSet(2:2)) == $TypeURSS{$K}((3:3, 5:6)))
                @test (setdiff(rs, 2:2) == $TypeURSS{$K}((3:3, 5:6)))
                @test (setdiff(rs, 2:2, 3:5) == $TypeURSS{$K}((6:6,)))

                rs1 = $TypeURSS{$K}((2:3, 5:6))
                rs2 = $TypeURSS{$K}((3:3, 5:6, 8:8))
                @inferred symdiff(rs1, rs2)
                @test collect(Iterators.flatten(symdiff(rs1, rs2))) ==
                      symdiff(collect(Iterators.flatten(rs1)), collect(Iterators.flatten(rs2)))
                @test collect(Iterators.flatten(symdiff(rs1, rs2, (9:9,)))) ==
                      symdiff(collect(Iterators.flatten(rs1)), collect(Iterators.flatten(rs2)), collect($K, Iterators.flatten((9:9,))))
                @test collect(Iterators.flatten(symdiff(rs1, (1:9,)))) ==
                      symdiff(collect(Iterators.flatten(rs1)), collect($K, Iterators.flatten((1:9,))))


                rs = $TypeURSS{$K}((0:0, 2:4))
                @test (union!(rs, $TypeURSS{$K}((0:0, 6:6))); rs == $TypeURSS{$K}((0:0, 2:4, 6:6)))
                @test pop!(rs, 0) == convertinfer($K, 0:0)
                @test_throws KeyError pop!(rs, 0)
                @test pop!(rs, 0, UInt64(1)) === UInt64(1)
                @test pop!(rs, 2, UInt64(1)) === convertinfer($K, 2:2)
                @test pop!(rs, 2, UInt64(1)) === UInt64(1)
                @test pop!(rs, 3:3) == convertinfer($K, 3:3)
                @test_throws KeyError pop!(rs, 3:3)
                @test pop!(rs, 3:3, UInt64(1)) === UInt64(1)
                @test rs == $TypeURSS{$K}((4:4, 6:6))
                @test pop!(rs, 4:4, UInt64(1)) == convertinfer($K, 4:4)
                @test pop!(rs, 4:4, UInt64(1)) === UInt64(1)
                @test size(rs) == (1,)
                @test isempty(rs) == false
                @test empty(rs) == $TypeURSS{$K}()
                empty!(rs)
                @test rs == $TypeURSS{$K}()
                @test size(rs) == (0,)
                @test isempty(rs) == true

                rs = $TypeURSS{$K}((0:0, 2:4))
                @test issubset(0, rs)
                @test issubset(0:0, rs)
                @test !issubset(0:1, rs)
                @test !issubset([0:1, 2:4], rs)
                @test issubset(rs, convertinfer($K, [0:1, 2:4]))
                @test issubset([0:0, 2:4], rs)
                @test issubset(rs, convertinfer($K, [0:0, 2:4]))
                @test !issubset(rs, convertinfer($K, [0:0, 3:4]))
                @test issubset(rs, convertinfer($K, 0:4))
                @test !issubset(rs, convertinfer($K, 0:3))
                rs2 = copy(rs)
                @test rs2 == rs && rs2 !== rs
                @test issubset(rs2, rs)
                @test issubset(rs, rs2)
                delete!(rs2, 0)
                @test issubset(rs2, rs)
                @test !issubset(rs, rs2)

            end
        end
    end
end


@testset "Iter allocations" begin
        for TypeURSS in list_of_containers_types_to_test

            @eval rs = $TypeURSS{UInt16}((1:6, 8:16))

            @test (@ballocated (for r in $rs j += 1 end) setup = (j = 0)) == 0
            @test (@ballocated (for r in $rs, i in r j += Int(i) end) setup = (j = 0)) == 0
            @test (@ballocated (for i in Iterators.flatten($rs) j += Int(i) end) setup = (j = 0)) == 0
            ss = subset(rs, 4:10)
            @test (@ballocated (for r in $ss j += 1 end) setup = (j = 0)) == 0
            @test (@ballocated (for r in $ss, i in r j += Int(i) end) setup = (j = 0)) == 0
            @test (@ballocated (for i in Iterators.flatten($ss) j += Int(i) end) setup = (j = 0)) == 0
            ss1 = subset(rs, 2:6)
            @test (@ballocated (for r in $ss1 j += 1 end) setup = (j = 0)) == 0
            @test (@ballocated (for r in $ss1, i in r j += Int(i) end) setup = (j = 0)) == 0
            @test (@ballocated (for i in Iterators.flatten($ss1) j += Int(i) end) setup = (j = 0)) == 0
            ss0 = subset(rs, 7:7)
            @test (@ballocated (for r in $ss0 j += 1 end) setup = (j = 0)) == 0
            @test (@ballocated (for r in $ss0, i in r j += Int(i) end) setup = (j = 0)) == 0
            @test (@ballocated (for i in Iterators.flatten($ss0) j += Int(i) end) setup = (j = 0)) == 0
            ss0 = subset(rs, 5:4)
            @test (@ballocated (for r in $ss0 j += 1 end) setup = (j = 0)) == 0
            @test (@ballocated (for r in $ss0, i in r j += Int(i) end) setup = (j = 0)) == 0
            @test (@ballocated (for i in Iterators.flatten($ss0) j += Int(i) end) setup = (j = 0)) == 0

        end
end


@testset "`show`" begin
    for K in list_of_Ti_to_test
        for TypeURSS in list_of_containers_types_to_test
            @eval begin

                io = IOBuffer()
                rs = $TypeURSS{$K}((0:0, 2:4))
                show(io, rs)
                TT = $TypeURSS
                KK = $K
                @test String(take!(io)) == "$TT{$KK}(" *
                                           repr($K(0)) * ":" * repr($K(0)) * ", " *
                                           repr($K(2)) * ":" * repr($K(4)) * ")"

                io = IOBuffer()
                rs = $TypeURSS{$K}((0:0, 2:4))
                show(io, MIME("text/plain"), rs)
                TT = $TypeURSS
                KK = $K
                pad = UnitRangesSortedSets.get_max_pad(rs)
                if length(rs) == 0
                    @test String(take!(io)) == "$TT{$KK}()\n"
                elseif length(rs) == 1
                    @test String(take!(io)) == "$TT{$KK} with 1 element:\n" *
                                               "  " * " "^(pad-length(repr($K(0)))) * repr($K(0)) * ":" * repr($K(0)) * "\n" *
                                               "  " * " "^(pad-length(repr($K(2)))) * repr($K(2)) * ":" * repr($K(4))
                else
                    @test String(take!(io)) == "$TT{$KK} with $(length(rs)) elements:\n" *
                                               "  " * " "^(pad-length(repr($K(0)))) * repr($K(0)) * ":" * repr($K(0)) * "\n" *
                                               "  " * " "^(pad-length(repr($K(2)))) * repr($K(2)) * ":" * repr($K(4))
                end

            end
        end
    end
end

@testset "Bugs" begin
    @testset "#5" begin
        for TypeURSS in list_of_containers_types_to_test
            @eval begin
                # provided by https://github.com/mnemnion
                tamil = Any['ꣳ', '௦':'௯', '௰':'௲', '௳', '𑌁', '𑌃', '𑌻':'𑌼', '𑿐':'𑿑', '𑿓', '᳚', '॒', '॑', '।', '॥']
                URvec = $TypeURSS[]
                for elem in tamil
                    push!(URvec, $TypeURSS([elem]))
                end
                @test_nowarn reduce(union, URvec)
                @test reduce(union, URvec) == $TypeURSS{Char}(('॑':'॒', '।':'॥', '௦':'௳', '᳚':'᳚', 'ꣳ':'ꣳ', '𑌁':'𑌁', '𑌃':'𑌃', '𑌻':'𑌼', '𑿐':'𑿑', '𑿓':'𑿓',))
            end
        end
    end
end
