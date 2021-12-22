using .UnitRangesSortedSets
using Test

# https://discourse.julialang.org/t/deparametrising-types/41939/4
basetype(t::DataType) = t.name.wrapper
basetype(t::UnionAll) = basetype(t.body)

function test_iseqial(rs1, rs2)
    length(rs1) == length(rs2) || return false
    for (r1, r2) in zip(rs1, rs2)
        r1 == r2 || return false
    end
    return true
end


@testset "Creating `Real`" begin
    for Ti in (Int, UInt64, UInt16, Float64)
        for TypeURSS in (UnitRangesSortedVector, UnitRangesSortedSet)
            @eval begin

                rs = $TypeURSS{$Ti}()
                @test eltype(rs) == UnitRange{$Ti}
                @test typeof(rs) === $TypeURSS{$Ti,UnitRange{$Ti}}
                @test length(rs) == 0
                @test test_iseqial(rs, Vector{UnitRange{$Ti}}(undef, 0))

                v = [1, 2]
                vu = [1:2]
                rs = $TypeURSS{$Ti}(v)
                @test eltype(rs) == UnitRange{$Ti}
                @test typeof(rs) === $TypeURSS{$Ti,UnitRange{$Ti}}
                @test length(rs) == 1
                @test test_iseqial(rs, vu)

                v = [1, 3, 4]
                vu = [1:1, 3:4]
                rs = $TypeURSS{$Ti}(v)
                @test eltype(rs) == UnitRange{$Ti}
                @test typeof(rs) === $TypeURSS{$Ti,UnitRange{$Ti}}
                @test length(rs) == 2
                @test test_iseqial(rs, vu)

                v = [$Ti(1), $Ti(3), $Ti(4)]
                vu = [1:1, 3:4]
                rs = $TypeURSS(v)
                @test eltype(rs) == UnitRange{$Ti}
                @test typeof(rs) === $TypeURSS{$Ti,UnitRange{$Ti}}
                @test length(rs) == 2
                @test test_iseqial(rs, vu)

                v = [3:4, 1:1]
                vu = [1:1, 3:4]
                rs = $TypeURSS{$Ti}(v)
                @test eltype(rs) == UnitRange{$Ti}
                @test typeof(rs) === $TypeURSS{$Ti,UnitRange{$Ti}}
                @test length(rs) == 2
                @test test_iseqial(rs, vu)

                v = [$Ti(3):$Ti(4), $Ti(1):$Ti(1)]
                vu = [1:1, 3:4]
                rs = $TypeURSS(v)
                @test eltype(rs) == UnitRange{$Ti}
                @test typeof(rs) === $TypeURSS{$Ti,UnitRange{$Ti}}
                @test length(rs) == 2
                @test test_iseqial(rs, vu)

                v = [3:4, 2:5, 1:1]
                vu = [1:5]
                rs = $TypeURSS{$Ti}(v)
                @test eltype(rs) == UnitRange{$Ti}
                @test typeof(rs) === $TypeURSS{$Ti,UnitRange{$Ti}}
                @test length(rs) == 1
                @test test_iseqial(rs, vu)

                v = [UInt64(3):UInt64(4), UInt64(1):UInt64(1)]
                vu = [UInt64(1):UInt64(1), UInt64(3):UInt64(4)]
                rs = $TypeURSS{$Ti}(v)
                @test eltype(rs) == UnitRange{$Ti}
                @test typeof(rs) === $TypeURSS{$Ti,UnitRange{$Ti}}
                @test length(rs) == 2
                @test test_iseqial(rs, vu)

                v = [UInt64(3), UInt64(4), UInt64(1)]
                vu = [UInt64(1):UInt64(1), UInt64(3):UInt64(4)]
                rs = $TypeURSS(v)
                @test eltype(rs) == UnitRange{UInt64}
                @test typeof(rs) === $TypeURSS{UInt64,UnitRange{UInt64}}
                @test length(rs) == 2
                @test test_iseqial(rs, vu)

                v = [UInt64(3):UInt64(4), UInt64(1):UInt64(1)]
                vu = [UInt64(1):UInt64(1), UInt64(3):UInt64(4)]
                rs = $TypeURSS(v)
                @test eltype(rs) == UnitRange{UInt64}
                @test typeof(rs) === $TypeURSS{UInt64,UnitRange{UInt64}}
                @test length(rs) == 2
                @test test_iseqial(rs, vu)
            end
        end
    end
end

@testset "Searching" begin
    for Ti in (Int, UInt64, UInt16, Float64, Char)
        for TypeURSS in (UnitRangesSortedVector, UnitRangesSortedSet)
            @eval begin
                rs = $TypeURSS{$Ti}((1:2, 4:4, 6:6))

                i = searchsortedlastrange(rs, 0)
                @test i == beforestartindex(rs)

                i = searchsortedlastrange(rs, 1)
                @test i != beforestartindex(rs) &&
                      getindex(rs, i) == UnitRange{$Ti}(1, 2) &&
                      i == firstindex(rs)

                i = searchsortedlastrange(rs, 2)
                @test i != beforestartindex(rs) &&
                      getindex(rs, i) == UnitRange{$Ti}(1, 2) &&
                      i == firstindex(rs)

                i = searchsortedlastrange(rs, 3)
                @test i != beforestartindex(rs) &&
                      getindex(rs, i) == UnitRange{$Ti}(1, 2) &&
                      i == firstindex(rs)

                i = searchsortedlastrange(rs, 4)
                @test i != beforestartindex(rs) &&
                      getindex(rs, i) == UnitRange{$Ti}(4, 4) &&
                      i != firstindex(rs) && i != lastindex(rs)

                i = searchsortedlastrange(rs, 5)
                @test i != beforestartindex(rs) &&
                      getindex(rs, i) == UnitRange{$Ti}(4, 4) &&
                      i != firstindex(rs) && i != lastindex(rs)

                i = searchsortedlastrange(rs, 6)
                @test i != beforestartindex(rs) &&
                      getindex(rs, i) == UnitRange{$Ti}(6, 6) &&
                      i == lastindex(rs)

                i = searchsortedlastrange(rs, 7)
                @test i != beforestartindex(rs) &&
                      getindex(rs, i) == UnitRange{$Ti}(6, 6) &&
                      i == lastindex(rs)


                i = searchsortedfirstrange(rs, 0)
                @test i != pastendindex(rs) &&
                      getindex(rs, i) == UnitRange{$Ti}(1, 2) &&
                      i == firstindex(rs)

                i = searchsortedfirstrange(rs, 1)
                @test i != pastendindex(rs) &&
                      getindex(rs, i) == UnitRange{$Ti}(1, 2) &&
                      i == firstindex(rs)

                i = searchsortedfirstrange(rs, 2)
                @test i != pastendindex(rs) &&
                      getindex(rs, i) == UnitRange{$Ti}(1, 2) &&
                      i == firstindex(rs)

                i = searchsortedfirstrange(rs, 3)
                @test i != pastendindex(rs) &&
                      getindex(rs, i) == UnitRange{$Ti}(4, 4) &&
                      i != firstindex(rs) && i != lastindex(rs)

                i = searchsortedfirstrange(rs, 4)
                @test i != pastendindex(rs) &&
                      getindex(rs, i) == UnitRange{$Ti}(4, 4) &&
                      i != firstindex(rs) && i != lastindex(rs)

                i = searchsortedfirstrange(rs, 5)
                @test i != pastendindex(rs) &&
                      getindex(rs, i) == UnitRange{$Ti}(6, 6) &&
                      i == lastindex(rs)

                i = searchsortedfirstrange(rs, 6)
                @test i != pastendindex(rs) &&
                      getindex(rs, i) == UnitRange{$Ti}(6, 6) &&
                      i == lastindex(rs)

                i = searchsortedfirstrange(rs, 7)
                @test i == pastendindex(rs)



                rs = $TypeURSS{$Ti}((1:2, 4:4))

                ii = searchsortedrange(rs, 0)
                @test length(ii) == 0 && getindex(rs, first(ii)) == UnitRange{$Ti}(1, 2) &&
                                         first(ii) == firstindex(rs) &&
                                         last(ii) == beforestartindex(rs)

                ii = searchsortedrange(rs, 1)
                @test length(ii) == 1 && getindex(rs, first(ii)) == UnitRange{$Ti}(1, 2) &&
                                         first(ii) == firstindex(rs)

                ii = searchsortedrange(rs, 2)
                @test length(ii) == 1 && getindex(rs, first(ii)) == UnitRange{$Ti}(1, 2)

                ii = searchsortedrange(rs, 3)
                @test length(ii) == 0 && getindex(rs, first(ii)) == UnitRange{$Ti}(4, 4) &&
                                         getindex(rs, last(ii)) == UnitRange{$Ti}(1, 2)

                ii = searchsortedrange(rs, 4)
                @test length(ii) == 1 && getindex(rs, first(ii)) == UnitRange{$Ti}(4, 4)

                ii = searchsortedrange(rs, 5)
                @test length(ii) == 0 && getindex(rs, last(ii)) == UnitRange{$Ti}(4, 4) &&
                                         first(ii) == pastendindex(rs) &&
                                         last(ii) == lastindex(rs)

                rs = $TypeURSS{$Ti}((1:2, 4:4))
                @test getrange(rs, 0) === nothing
                @test getrange(rs, 1) == UnitRange{$Ti}(1, 1)
                @test getrange(rs, 1:1) == UnitRange{$Ti}(1, 1)
                @test getrange(rs, 1:2) == UnitRange{$Ti}(1, 2)
                @test getrange(rs, 1:3) === nothing
            end
        end
    end
end

@testset "`in`" begin
    for Ti in (Int, UInt64, UInt16, Float64, Char)
        for TypeURSS in (UnitRangesSortedVector, UnitRangesSortedSet)
            @eval begin

                v = [1, 2, 3, 5, 6, 8]
                rs = $TypeURSS{$Ti}(v)
                for i = 0:9
                    @test !xor(in(i, rs), in(i, v))
                end

                v = [2:3, 5:6]
                rs = $TypeURSS{$Ti}(v)
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
    for Ti in (Int, UInt64, UInt16, Float64)
        for TypeURSS in (UnitRangesSortedVector, UnitRangesSortedSet)
            @eval begin

                rs = $TypeURSS{$Ti}()
                @test test_iseqial(push!(rs, 3),  (3:3,))
                @test test_iseqial(push!(rs, 4),  (3:4,))
                @test test_iseqial(push!(rs, 18), (3:4, 18:18))
                @test test_iseqial(push!(rs, 2),  (2:4, 18:18))
                @test test_iseqial(push!(rs, 0),  (0:0, 2:4, 18:18))
                @test test_iseqial(push!(rs, 20), (0:0, 2:4, 18:18, 20:20))
                @test test_iseqial(push!(rs, 19), (0:0, 2:4, 18:20))
                @test test_iseqial(push!(rs, 1),  (0:4, 18:20))
                @test test_iseqial(push!(rs, 9),  (0:4, 9:9, 18:20))
                @test test_iseqial(push!(rs, 9),  (0:4, 9:9, 18:20))
                @test test_iseqial(push!(rs, 10), (0:4, 9:10, 18:20))
                @test test_iseqial(push!(rs, 8),  (0:4, 8:10, 18:20))

                rs = $TypeURSS{$Ti}()
                @test test_iseqial(push!(rs, 13:14), (13:14,))
                @test test_iseqial(push!(rs, 23:24), (13:14, 23:24))
                @test test_iseqial(push!(rs, 12:12), (12:14, 23:24))
                @test test_iseqial(push!(rs, 11:13), (11:14, 23:24))
                @test test_iseqial(push!(rs, 8:9),   (8:9, 11:14, 23:24))
                @test test_iseqial(push!(rs, 9:12),  (8:14, 23:24))
                @test test_iseqial(push!(rs, 16:17),  (8:14, 16:17, 23:24))
                @test test_iseqial(push!(rs, 25:26),  (8:14, 16:17, 23:26))
                @test test_iseqial(push!(rs, 24:26),  (8:14, 16:17, 23:26))
                @test test_iseqial(push!(rs, 28:29),  (8:14, 16:17, 23:26, 28:29))
                @test test_iseqial(push!(rs, 28:29),  (8:14, 16:17, 23:26, 28:29))
                @test test_iseqial(push!(rs, 16:15),  (8:14, 16:17, 23:26, 28:29))
                @test test_iseqial(push!(rs, 18:19),  (8:14, 16:19, 23:26, 28:29))
            end
        end
    end
end


@testset "`delete!`" begin
    for Ti in (Int, UInt64, UInt16, Float64)
        for TypeURSS in (UnitRangesSortedVector, UnitRangesSortedSet)
            @eval begin

                rs = $TypeURSS{$Ti}((0:0, 2:4, 8:8, 10:12))
                @test test_iseqial(delete!(rs, 5),  (0:0, 2:4, 8:8, 10:12))
                @test test_iseqial(delete!(rs, 0),  (2:4, 8:8, 10:12))
                @test test_iseqial(delete!(rs, 0),  (2:4, 8:8, 10:12))
                @test test_iseqial(delete!(rs, 8),  (2:4, 10:12))
                @test test_iseqial(delete!(rs, 12), (2:4, 10:11))
                @test test_iseqial(delete!(rs, 10), (2:4, 11:11))
                @test test_iseqial(delete!(rs, 11), (2:4,))
                @test test_iseqial(delete!(rs, 3),  (2:2, 4:4))
                @test delete!(rs, 3) === rs

                rs = $TypeURSS{$Ti}((0:0, 2:6, 8:8, 10:13, 15:20))
                @test test_iseqial(delete!(rs, 17:17),  (0:0, 2:6, 8:8, 10:13, 15:16, 18:20))
                @test test_iseqial(delete!(rs, 16:18),  (0:0, 2:6, 8:8, 10:13, 15:15, 19:20))
                @test test_iseqial(delete!(rs, 15:20),  (0:0, 2:6, 8:8, 10:13))
                @test test_iseqial(delete!(rs, 7:7),  (0:0, 2:6, 8:8, 10:13))
                @test test_iseqial(delete!(rs, 6:7),  (0:0, 2:5, 8:8, 10:13))
                @test test_iseqial(delete!(rs, 0:1),  (2:5, 8:8, 10:13))
                @test test_iseqial(delete!(rs, 9:10),  (2:5, 8:8, 11:13))
                @test test_iseqial(delete!(rs, 13:14),  (2:5, 8:8, 11:12))
                @test test_iseqial(delete!(rs, 1:3),  (4:5, 8:8, 11:12))

            end
        end
    end
end


function test_iterators(rs, vu, T::Type)
    @test basetype(typeof(rs)) == T
    @test length(rs) == length(vu)
    @test rs[begin] == first(rs) == first(vu)
    @test rs[end] == last(rs) == last(vu)
    @test test_iseqial(rs, vu)
    @test test_iseqial(Iterators.reverse(rs), Iterators.reverse(vu))
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
    @test test_iseqial(rs, vu)
    @test test_iseqial(Iterators.reverse(rs), Iterators.reverse(vu))
    @test [r for r in rs] == [r for r in vu]
    @test [r for r in Iterators.reverse(rs)] == [r for r in Iterators.reverse(vu)]
    @test [x for r in rs for x in r] == [x for r in vu for x in r]
    @test [x for r in Iterators.reverse(rs) for x in r] == [x for r in Iterators.reverse(vu) for x in r]
    @test_throws ArgumentError sum(r->length(r), rs)
end


@testset "Common functions" begin

    for Ti in (Int, UInt64, UInt16, Float64)
        for TypeURSS in (UnitRangesSortedVector, UnitRangesSortedSet)
            @eval begin

                rs = $TypeURSS{$Ti}()
                vu = Vector{eltype(ss)}([])
                test_empty_iterators(rs, vu, $TypeURSS)

                rs = $TypeURSS{$Ti}((1:6,))
                vu = (1:6,)
                test_iterators(rs, vu, $TypeURSS)

                rs = $TypeURSS{$Ti}((1:6, 8:16, 20:33))
                vu = (1:6, 8:16, 20:33)
                test_iterators(rs, vu, $TypeURSS)

                rs = $TypeURSS{$Ti}((0:0, 2:4))
                rs2 = copy(rs)
                @test rs2 == rs && rs2 !== rs
                rs2 = empty(rs)
                @test rs2 == $TypeURSS{$Ti}()
                empty!(rs2)
                @test rs2 != rs && rs2 !== rs
                @test length(rs2) == 0

                rs = $TypeURSS{$Ti}((0:0, 2:4))
                @test filter(s->length(s)==3, rs) == $TypeURSS{$Ti}((2:4))
                @test (filter!(s->length(s)==1, rs); rs == $TypeURSS{$Ti}((0:0)))

                rs = $TypeURSS{$Ti}((0:0, 2:4))
                @test collect(rs) == Vector{basetype(eltype(rs)){UInt64}}([0:0, 2:4])
                @test collect(UInt64, rs) == Vector{basetype(eltype(rs)){UInt64}}([0:0, 2:4])
                @test convert(typeof(rs), rs) === rs
                @test convert(Vector{$Ti}, rs) == Vector{$Ti}([0,2,3,4])
                @test convert(Set{$Ti}, rs) == Set{$Ti}([0,2,3,4])

            end
        end
    end
end


@testset "`subset`" begin

    for Ti in (Int, UInt64, UInt16, Float64)
        for TypeURSS in (UnitRangesSortedVector, UnitRangesSortedSet)
            @eval begin

                rs = $TypeURSS{$Ti}((1:6, 8:16, 20:33, 35:47, 49:50))

                ss = subset(rs, 1:50)
                vu = (1:6, 8:16, 20:33, 35:47, 49:50)
                test_iterators(ss, vu, SubUnitRangesSortedSet)

                ss = subset(rs, 10:40)
                vu = (10:16, 20:33, 35:40)
                test_iterators(ss, vu, SubUnitRangesSortedSet)

                ss = subset(rs, 7:40)
                vu = (8:16, 20:33, 35:40)
                test_iterators(ss, vu, SubUnitRangesSortedSet)

                ss = subset(rs, 11:48)
                vu = (11:16, 20:33, 35:47)
                test_iterators(ss, vu, SubUnitRangesSortedSet)

                ss = subset(rs, 7:48)
                vu = (8:16, 20:33, 35:47)
                test_iterators(ss, vu, SubUnitRangesSortedSet)

                ss = subset(rs, 21:32)
                vu = (21:32,)
                test_iterators(ss, vu, Sub1UnitRangesSortedSet)

                ss = subset(rs, 20:33)
                vu = (20:33,)
                test_iterators(ss, vu, Sub1UnitRangesSortedSet)

                ss = subset(rs, 19:34)
                vu = (20:33,)
                test_iterators(ss, vu, Sub1UnitRangesSortedSet)

                ss = subset(rs, 22:34)
                vu = (22:33,)
                test_iterators(ss, vu, Sub1UnitRangesSortedSet)

                ss = subset(rs, 19:30)
                vu = (20:30,)
                test_iterators(ss, vu, Sub1UnitRangesSortedSet)

                ss = subset(rs, 18:19)
                vu = Vector{eltype(ss)}([])
                test_empty_iterators(ss, vu, Sub0UnitRangesSortedSet)

                ss = subset(rs, 21:20)
                vu = Vector{eltype(ss)}([])
                test_empty_iterators(ss, vu, Sub0UnitRangesSortedSet)

            end
        end
    end
end


@testset "`Set` API" begin
    for Ti in (Int, UInt64, UInt16, Float64)
        for TypeURSS in (UnitRangesSortedVector, UnitRangesSortedSet)
            @eval begin

                @test union($TypeURSS{$Ti}((0:0, 2:4)), $TypeURSS{$Ti}((2:3, 5:6))) == $TypeURSS{$Ti}((0:0, 2:6))
                @test union($TypeURSS{$Ti}((0:0, 2:4)), (2:3, 5:6)) == $TypeURSS{$Ti}((0:0, 2:6))
                @test union($TypeURSS{$Ti}((0:0, 2:4)), (2, 3, 5, 6)) == $TypeURSS{$Ti}((0:0, 2:6))
                @test union($TypeURSS{$Ti}((0:0, 2:4)), [2:3, 5:6]) == $TypeURSS{$Ti}((0:0, 2:6))
                @test union($TypeURSS{$Ti}((0:0, 2:4)), [2, 3, 5, 6]) == $TypeURSS{$Ti}((0:0, 2:6))
                @test union($TypeURSS{$Ti}((0:0, 2:4)), 1) == $TypeURSS{$Ti}((0:4))

                @test intersect($TypeURSS{$Ti}((0:0, 2:4)), $TypeURSS{$Ti}((2:3, 5:6))) == $TypeURSS{$Ti}((2:3))
                @test intersect($TypeURSS{$Ti}((0:0, 2:4)), (2:3, 5:6)) == $TypeURSS{$Ti}((2:3))
                @test intersect($TypeURSS{$Ti}((0:0, 2:4)), 2) == $TypeURSS{$Ti}((2:2))
                rs = $TypeURSS{$Ti}((0:0, 2:4))
                @test (intersect!(rs, $TypeURSS{$Ti}((2:3, 5:6, 8:8))); rs == $TypeURSS{$Ti}((2:3)))
                rs = [0:0, 2:4]
                @test (intersect!(rs, $TypeURSS{$Ti}((2:3, 5:6, 8:8))); rs == [2:3])

                rs = $TypeURSS{$Ti}((0:0, 2:4))
                @test (union!(rs, $TypeURSS{$Ti}((0:0, 6:6))); rs == $TypeURSS{$Ti}((0:0, 2:4, 6:6)))
                @test pop!(rs, 0) == UnitRange{$Ti}(0, 0)
                @test_throws KeyError pop!(rs, 0)
                @test pop!(rs, 0, UInt64(1)) === UInt64(1)
                @test pop!(rs, 2, UInt64(1)) === UnitRange{$Ti}(2, 2)
                @test pop!(rs, 2, UInt64(1)) === UInt64(1)
                @test pop!(rs, 3:3) == UnitRange{$Ti}(3, 3)
                @test_throws KeyError pop!(rs, 3:3)
                @test pop!(rs, 3:3, UInt64(1)) === UInt64(1)
                @test rs == $TypeURSS{$Ti}((4:4, 6:6))
                @test pop!(rs, 4:4, UInt64(1)) == UnitRange{$Ti}(4, 4)
                @test pop!(rs, 4:4, UInt64(1)) === UInt64(1)
                @test empty(rs) == $TypeURSS{$Ti}()
                empty!(rs)
                @test rs == $TypeURSS{$Ti}()

                rs = $TypeURSS{$Ti}((0:0, 2:4))
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
    for Ti in (Int, UInt64, UInt16, Float64)
        for TypeURSS in (UnitRangesSortedVector, UnitRangesSortedSet)
            @eval begin

                io = IOBuffer()
                rs = $TypeURSS{$Ti}((0:0, 2:4))
                show(io, rs)
                T = $TypeURSS{$Ti}
                @test String(take!(io)) == "$T(" *
                                           repr($Ti(0)) * ":" * repr($Ti(0)) * ", " *
                                           repr($Ti(2)) * ":" * repr($Ti(4)) * ")"

                io = IOBuffer()
                rs = $TypeURSS{$Ti}((0:0, 2:4))
                show(io, MIME("text/plain"), rs)
                T = $TypeURSS{$Ti}
                @test String(take!(io)) == "$T():\n" *
                                           "  " * repr($Ti(0)) * ":" * repr($Ti(0)) * "\n" *
                                           "  " * repr($Ti(2)) * ":" * repr($Ti(4))

            end
        end
    end
end


