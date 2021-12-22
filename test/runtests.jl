using .UnitRangesSortedSets
import .UnitRangesSortedSets: inferrangetype, to_urange
using Test

# https://github.com/JuliaLang/julia/issues/39952
basetype(::Type{T}) where T = Base.typename(T).wrapper


const list_of_Ti_to_test = (Int, UInt64, UInt16, Float64, Char)

const list_of_containers_types_to_test = (UnitRangesSortedVector, UnitRangesSortedSet)


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
    Trange = inferrangetype(Tv)
    @test eltype(rs) === Trange
    @test typeof(rs) === basetype(T){Tv,Trange}
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

                ir = searchsortedlastrange(rs, 0)
                @test ir == beforestartindex(rs)

                ir = searchsortedlastrange(rs, 1)
                @test ir != beforestartindex(rs) &&
                      getindex(rs, ir) == convertinfer($K, 1:2) &&
                      ir == firstindex(rs)

                ir = searchsortedlastrange(rs, 2)
                @test ir != beforestartindex(rs) &&
                      getindex(rs, ir) == convertinfer($K, 1:2) &&
                      ir == firstindex(rs)

                ir = searchsortedlastrange(rs, 3)
                @test ir != beforestartindex(rs) &&
                      getindex(rs, ir) == convertinfer($K, 1:2) &&
                      ir == firstindex(rs)

                ir = searchsortedlastrange(rs, 4)
                @test ir != beforestartindex(rs) &&
                      getindex(rs, ir) == convertinfer($K, 4:4) &&
                      ir != firstindex(rs) && ir != lastindex(rs)

                ir = searchsortedlastrange(rs, 5)
                @test ir != beforestartindex(rs) &&
                      getindex(rs, ir) == convertinfer($K, 4:4) &&
                      ir != firstindex(rs) && ir != lastindex(rs)

                ir = searchsortedlastrange(rs, 6)
                @test ir != beforestartindex(rs) &&
                      getindex(rs, ir) == convertinfer($K, 6:6) &&
                      ir == lastindex(rs)

                ir = searchsortedlastrange(rs, 7)
                @test ir != beforestartindex(rs) &&
                      getindex(rs, ir) == convertinfer($K, 6:6) &&
                      ir == lastindex(rs)


                ir = searchsortedfirstrange(rs, 0)
                @test ir != pastendindex(rs) &&
                      getindex(rs, ir) == convertinfer($K, 1:2) &&
                      ir == firstindex(rs)

                ir = searchsortedfirstrange(rs, 1)
                @test ir != pastendindex(rs) &&
                      getindex(rs, ir) == convertinfer($K, 1:2) &&
                      ir == firstindex(rs)

                ir = searchsortedfirstrange(rs, 2)
                @test ir != pastendindex(rs) &&
                      getindex(rs, ir) == convertinfer($K, 1:2) &&
                      ir == firstindex(rs)

                ir = searchsortedfirstrange(rs, 3)
                @test ir != pastendindex(rs) &&
                      getindex(rs, ir) == convertinfer($K, 4:4) &&
                      ir != firstindex(rs) && ir != lastindex(rs)

                ir = searchsortedfirstrange(rs, 4)
                @test ir != pastendindex(rs) &&
                      getindex(rs, ir) == convertinfer($K, 4:4) &&
                      ir != firstindex(rs) && ir != lastindex(rs)

                ir = searchsortedfirstrange(rs, 5)
                @test ir != pastendindex(rs) &&
                      getindex(rs, ir) == convertinfer($K, 6:6) &&
                      ir == lastindex(rs)

                ir = searchsortedfirstrange(rs, 6)
                @test ir != pastendindex(rs) &&
                      getindex(rs, ir) == convertinfer($K, 6:6) &&
                      ir == lastindex(rs)

                ir = searchsortedfirstrange(rs, 7)
                @test ir == pastendindex(rs)



                rs = $TypeURSS{$K}((1:2, 4:4))

                ii = searchsortedrange(rs, 0)
                @test length(ii) == 0 && getindex(rs, first(ii)) == convertinfer($K, 1:2) &&
                                         first(ii) == firstindex(rs) &&
                                         last(ii) == beforestartindex(rs)

                ii = searchsortedrange(rs, 1)
                @test length(ii) == 1 && getindex(rs, first(ii)) == convertinfer($K, 1:2) &&
                                         first(ii) == firstindex(rs)

                ii = searchsortedrange(rs, 2)
                @test length(ii) == 1 && getindex(rs, first(ii)) == convertinfer($K, 1:2)

                ii = searchsortedrange(rs, 3)
                @test length(ii) == 0 && getindex(rs, first(ii)) == convertinfer($K, 4:4) &&
                                         getindex(rs, last(ii)) == convertinfer($K, 1:2)

                ii = searchsortedrange(rs, 4)
                @test length(ii) == 1 && getindex(rs, first(ii)) == convertinfer($K, 4:4)

                ii = searchsortedrange(rs, 5)
                @test length(ii) == 0 && getindex(rs, last(ii)) == convertinfer($K, 4:4) &&
                                         first(ii) == pastendindex(rs) &&
                                         last(ii) == lastindex(rs)

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
                @test check_isequal(delete!(rs, 1:3),  convertinfer($K, (4:5, 8:8, 11:12)))

            end
        end
    end
end


function test_iterators(rs, vu, T::Type)
    @test basetype(typeof(rs)) == T
    @test length(rs) == length(vu)
    @test check_ranges_isequal(rs[begin], first(rs), first(vu))
    @test check_ranges_isequal(rs[end], last(rs), last(vu))
    @test check_isequal(rs, vu)
    @test check_isequal(Iterators.reverse(rs), Iterators.reverse(vu))
    @test [r for r in rs] == [r for r in vu]
    @test [r for r in Iterators.reverse(rs)] == [r for r in Iterators.reverse(vu)]
    @test [x for r in rs for x in r] == [x for r in vu for x in r]
    @test [x for r in Iterators.reverse(rs) for x in r] == [x for r in Iterators.reverse(vu) for x in r]
    @test sum(r->length(r), rs) ==
          sum(r->length(r), Iterators.reverse(rs)) ==
          sum(i->length(rs[i]), eachindex(rs)) ==
          sum(i->length(rs[i]), Iterators.reverse(eachindex(rs))) ==
          sum(r->length(r), vu)
end

function test_empty_iterators(rs, vu, T::Type)
    @test basetype(typeof(rs)) == T
    @test length(rs) == length(vu)
    @test_throws BoundsError rs[begin]
    @test_throws BoundsError rs[end]
    @test check_isequal(rs, vu)
    @test check_isequal(Iterators.reverse(rs), Iterators.reverse(vu))
    @test [r for r in rs] == [r for r in vu]
    @test [r for r in Iterators.reverse(rs)] == [r for r in Iterators.reverse(vu)]
    @test [x for r in rs for x in r] == [x for r in vu for x in r]
    @test [x for r in Iterators.reverse(rs) for x in r] == [x for r in Iterators.reverse(vu) for x in r]
    @test_throws ArgumentError sum(r->length(r), rs)
end


@testset "Common functions" begin
    for K in list_of_Ti_to_test
        for TypeURSS in list_of_containers_types_to_test
            @eval begin

                rs = $TypeURSS{$K}()
                vu = Vector{eltype(rs)}([])
                test_empty_iterators(rs, convertinfer($K, vu), $TypeURSS)

                rs = $TypeURSS{$K}((1:6,))
                vu = (1:6,)
                test_iterators(rs, convertinfer($K, vu), $TypeURSS)

                rs = $TypeURSS{$K}((1:6, 8:16, 20:33))
                vu = (1:6, 8:16, 20:33)
                test_iterators(rs, convertinfer($K, vu), $TypeURSS)

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
                @test collect(rs) == convertinfer($K, [0:0, 2:4])
                @test collect(UInt64, rs) == convertinfer(UInt64, [0:0, 2:4])
                @test convert(typeof(rs), rs) === rs
                @test convert(Vector{$K}, rs) == Vector{$K}([0,2,3,4])
                @test convert(Vector{UInt64}, rs) == Vector{UInt64}([0,2,3,4])
                @test convert(Set{$K}, rs) == Set{$K}(convertinfer($K, [0,2,3,4]))
                @test convert(Set{UInt64}, rs) == Set{UInt64}(convertinfer(UInt64, [0,2,3,4]))

            end
        end
    end
end


@testset "`subset`" begin

    for K in (Int, UInt64, UInt16, Float64)
        for TypeURSS in (UnitRangesSortedVector, UnitRangesSortedSet)
            @eval begin

                rs = $TypeURSS{$K}((1:6, 8:16, 20:33, 35:47, 49:50))

                ss = subset(rs, 1:50)
                vu = (1:6, 8:16, 20:33, 35:47, 49:50)
                test_iterators(ss, convertinfer($K, vu), SubUnitRangesSortedSet)

                ss = subset(rs, 10:40)
                vu = (10:16, 20:33, 35:40)
                test_iterators(ss, convertinfer($K, vu), SubUnitRangesSortedSet)

                ss = subset(rs, 7:40)
                vu = (8:16, 20:33, 35:40)
                test_iterators(ss, convertinfer($K, vu), SubUnitRangesSortedSet)

                ss = subset(rs, 11:48)
                vu = (11:16, 20:33, 35:47)
                test_iterators(ss, convertinfer($K, vu), SubUnitRangesSortedSet)

                ss = subset(rs, 7:48)
                vu = (8:16, 20:33, 35:47)
                test_iterators(ss, convertinfer($K, vu), SubUnitRangesSortedSet)

                ss = subset(rs, 21:32)
                vu = (21:32,)
                test_iterators(ss, convertinfer($K, vu), Sub1UnitRangesSortedSet)

                ss = subset(rs, 20:33)
                vu = (20:33,)
                test_iterators(ss, convertinfer($K, vu), Sub1UnitRangesSortedSet)

                ss = subset(rs, 19:34)
                vu = (20:33,)
                test_iterators(ss, convertinfer($K, vu), Sub1UnitRangesSortedSet)

                ss = subset(rs, 22:34)
                vu = (22:33,)
                test_iterators(ss, convertinfer($K, vu), Sub1UnitRangesSortedSet)

                ss = subset(rs, 19:30)
                vu = (20:30,)
                test_iterators(ss, convertinfer($K, vu), Sub1UnitRangesSortedSet)

                ss = subset(rs, 18:19)
                vu = Vector{eltype(ss)}([])
                test_empty_iterators(ss, convertinfer($K, vu), Sub0UnitRangesSortedSet)

                ss = subset(rs, 21:20)
                vu = Vector{eltype(ss)}([])
                test_empty_iterators(ss, convertinfer($K, vu), Sub0UnitRangesSortedSet)

            end
        end
    end
end


@testset "`Set` API" begin
    for K in (Int, UInt64, UInt16, Float64)
        for TypeURSS in (UnitRangesSortedVector, UnitRangesSortedSet)
            @eval begin

                @test union($TypeURSS{$K}((0:0, 2:4)), $TypeURSS{$K}((2:3, 5:6))) == $TypeURSS{$K}((0:0, 2:6))
                @test union($TypeURSS{$K}((0:0, 2:4)), (2:3, 5:6)) == $TypeURSS{$K}((0:0, 2:6))
                @test union($TypeURSS{$K}((0:0, 2:4)), (2, 3, 5, 6)) == $TypeURSS{$K}((0:0, 2:6))
                @test union($TypeURSS{$K}((0:0, 2:4)), [2:3, 5:6]) == $TypeURSS{$K}((0:0, 2:6))
                @test union($TypeURSS{$K}((0:0, 2:4)), [2, 3, 5, 6]) == $TypeURSS{$K}((0:0, 2:6))
                @test union($TypeURSS{$K}((0:0, 2:4)), 1) == $TypeURSS{$K}((0:4))

                @test intersect($TypeURSS{$K}((0:0, 2:4)), $TypeURSS{$K}((2:3, 5:6))) == $TypeURSS{$K}((2:3))
                @test intersect($TypeURSS{$K}((0:0, 2:4)), (2:3, 5:6)) == $TypeURSS{$K}((2:3))
                @test intersect($TypeURSS{$K}((0:0, 2:4)), 2) == $TypeURSS{$K}((2:2))
                rs = $TypeURSS{$K}((0:0, 2:4))
                @test (intersect!(rs, $TypeURSS{$K}((2:3, 5:6, 8:8))); rs == $TypeURSS{$K}((2:3)))
                rs = [0:0, 2:4]
                @test (intersect!(rs, $TypeURSS{$K}((2:3, 5:6, 8:8))); rs == [2:3])

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
                @test empty(rs) == $TypeURSS{$K}()
                empty!(rs)
                @test rs == $TypeURSS{$K}()

                rs = $TypeURSS{$K}((0:0, 2:4))
                @test issubset(0, rs)
                @test issubset(0:0, rs)
                @test !issubset(0:1, rs)
                @test !issubset([0:1, 2:4], rs)
                @test issubset(rs, [0:1, 2:4])
                @test issubset([0:0, 2:4], rs)
                @test issubset(rs, [0:0, 2:4])
                @test !issubset(rs, [0:0, 3:4])
                @test issubset(rs, 0:4)
                @test !issubset(rs, 0:3)
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


@testset "`show`" begin
    for K in (Int, UInt64, UInt16, Float64)
        for TypeURSS in (UnitRangesSortedVector, UnitRangesSortedSet)
            @eval begin

                io = IOBuffer()
                rs = $TypeURSS{$K}((0:0, 2:4))
                show(io, rs)
                T = $TypeURSS{$K}
                @test String(take!(io)) == "$T(" *
                                           repr($K(0)) * ":" * repr($K(0)) * ", " *
                                           repr($K(2)) * ":" * repr($K(4)) * ")"

                io = IOBuffer()
                rs = $TypeURSS{$K}((0:0, 2:4))
                show(io, MIME("text/plain"), rs)
                T = $TypeURSS{$K}
                @test String(take!(io)) == "$T():\n" *
                                           "  " * repr($K(0)) * ":" * repr($K(0)) * "\n" *
                                           "  " * repr($K(2)) * ":" * repr($K(4))

            end
        end
    end
end


