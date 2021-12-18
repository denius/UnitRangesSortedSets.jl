using .UnitRangesSortedSets
using Test

function test_iseqial(rs1, rs2)
    length(rs1) == length(rs2) || return false
    for (r1, r2) in zip(rs1, rs2)
        r1 == r2 || return false
    end
    return true
end


@testset "creating" begin
    for Ti in (Int, UInt, UInt16, Float64)
        for TypeURSS in (UnitRangesSortedVector, UnitRangesSortedSet)
            @eval begin

                rs = $TypeURSS{$Ti}()
                @test eltype(rs) == $Ti
                @test typeof(rs) === $TypeURSS{$Ti}
                @test length(rs) == 0
                @test test_iseqial(rs, Vector{UnitRange{$Ti}}(undef, 0))

                v = [1, 2]
                vu = [1:2]
                rs = $TypeURSS{$Ti}(v)
                @test eltype(rs) == $Ti
                @test typeof(rs) === $TypeURSS{$Ti}
                @test length(rs) == 1
                @test test_iseqial(rs, vu)

                v = [1, 3, 4]
                vu = [1:1, 3:4]
                rs = $TypeURSS{$Ti}(v)
                @test eltype(rs) == $Ti
                @test typeof(rs) === $TypeURSS{$Ti}
                @test length(rs) == 2
                @test test_iseqial(rs, vu)

                v = [3:4, 1:1]
                vu = [1:1, 3:4]
                rs = $TypeURSS{$Ti}(v)
                @test eltype(rs) == $Ti
                @test typeof(rs) === $TypeURSS{$Ti}
                @test length(rs) == 2
                @test test_iseqial(rs, vu)

                v = [3:4, 2:5, 1:1]
                vu = [1:5]
                rs = $TypeURSS{$Ti}(v)
                @test eltype(rs) == $Ti
                @test typeof(rs) === $TypeURSS{$Ti}
                @test length(rs) == 1
                @test test_iseqial(rs, vu)

                v = [UInt(3):UInt(4), UInt(1):UInt(1)]
                vu = [UInt(1):UInt(1), UInt(3):UInt(4)]
                rs = $TypeURSS{$Ti}(v)
                @test eltype(rs) == $Ti
                @test typeof(rs) === $TypeURSS{$Ti}
                @test length(rs) == 2
                @test test_iseqial(rs, vu)

                v = [UInt(3):UInt(4), UInt(1):UInt(1)]
                vu = [UInt(1):UInt(1), UInt(3):UInt(4)]
                rs = $TypeURSS(v)
                @test eltype(rs) == UInt
                @test typeof(rs) === $TypeURSS{UInt}
                @test length(rs) == 2
                @test test_iseqial(rs, vu)
            end
        end
    end
end


@testset "in" begin
    for Ti in (Int, UInt, UInt16, Float64)
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
                @test !in(2:4, rs)
                @test !in(4:4, rs)
                @test !in(2:6, rs)
                @test !in(1:7, rs)
            end
        end
    end
end


@testset "push!" begin
    for Ti in (Int, UInt, UInt16, Float64)
        for TypeURSS in (UnitRangesSortedVector, UnitRangesSortedSet)
            @eval begin

                rs = $TypeURSS{$Ti}()
                @test test_iseqial(push!(rs, 3),  (3:3,))
                @test test_iseqial(push!(rs, 4),  (3:4,))
                @test test_iseqial(push!(rs, 8),  (3:4, 8:8))
                @test test_iseqial(push!(rs, 2),  (2:4, 8:8))
                @test test_iseqial(push!(rs, 0),  (0:0, 2:4, 8:8))
                @test test_iseqial(push!(rs, 10), (0:0, 2:4, 8:8, 10:10))
                @test test_iseqial(push!(rs, 9),  (0:0, 2:4, 8:10))
                @test test_iseqial(push!(rs, 1),  (0:4, 8:10))
                @test test_iseqial(push!(rs, 6),  (0:4, 6:6, 8:10))
                @test test_iseqial(push!(rs, 6),  (0:4, 6:6, 8:10))

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
            end
        end
    end
end


@testset "delete!" begin
    for Ti in (Int, UInt, UInt16, Float64)
        for TypeURSS in (UnitRangesSortedVector, UnitRangesSortedSet)
            @eval begin

                rs = $TypeURSS{$Ti}((0:0, 2:4, 8:8, 10:12))
                @test test_iseqial(delete!(rs, 5),  (0:0, 2:4, 8:8, 10:12))
                @test test_iseqial(delete!(rs, 0),  (2:4, 8:8, 10:12))
                @test test_iseqial(delete!(rs, 8),  (2:4, 10:12))
                @test test_iseqial(delete!(rs, 12), (2:4, 10:11))
                @test test_iseqial(delete!(rs, 10), (2:4, 11:11))
                @test test_iseqial(delete!(rs, 11), (2:4,))
                @test test_iseqial(delete!(rs, 3),  (2:2, 4:4))
                @test delete!(rs, 3) === rs

                rs = $TypeURSS{$Ti}((0:0, 2:6, 8:8, 10:13))
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


@testset "Some functions" begin
    for Ti in (Int, UInt, UInt16, Float64)
        for TypeURSS in (UnitRangesSortedVector, UnitRangesSortedSet)
            @eval begin

                @test union($TypeURSS{$Ti}((0:0, 2:4)), $TypeURSS{$Ti}((2:3, 5:6))) == $TypeURSS{$Ti}((0:0, 2:6))
                @test union($TypeURSS{$Ti}((0:0, 2:4)), (2:3, 5:6)) == $TypeURSS{$Ti}((0:0, 2:6))
                @test union($TypeURSS{$Ti}((0:0, 2:4)), 1) == $TypeURSS{$Ti}((0:4))

                @test intersect($TypeURSS{$Ti}((0:0, 2:4)), $TypeURSS{$Ti}((2:3, 5:6))) == $TypeURSS{$Ti}((2:3))
                @test intersect($TypeURSS{$Ti}((0:0, 2:4)), (2:3, 5:6)) == $TypeURSS{$Ti}((2:3))
                @test intersect($TypeURSS{$Ti}((0:0, 2:4)), 2) == $TypeURSS{$Ti}((2:2))
                rs = $TypeURSS{$Ti}((0:0, 2:4))
                @test (intersect!(rs, $TypeURSS{$Ti}((2:3, 5:6))); rs == $TypeURSS{$Ti}((2:3)))

                rs = $TypeURSS{$Ti}((0:0, 2:4))
                @test (union!(rs, $TypeURSS{$Ti}((0:0, 6:6))); rs == $TypeURSS{$Ti}((0:0, 2:4, 6:6)))
                @test pop!(rs, 0) == Ti(0)
                @test_throws KeyError pop!(rs, 0)
                @test pop!(rs, 0, Ti(0)) == Ti(0)
                empty!(rs)
                @test rs == $TypeURSS{$Ti}()

                rs = $TypeURSS{$Ti}((0:0, 2:4))
                @test issubset(0, rs)
                @test issubset(0:0, rs)
                @test !issubset(0:1, rs)
                @test !issubset([0:1, 2:4], rs)
                @test issubset([0:0, 2:4], rs)
                rs2 = copy(rs)
                @test rs2 == rs && rs2 !== rs
                rs2 = copymutable(rs)
                rs2 = emptymutable(rs)
                @test rs2 == $TypeURSS{$Ti}()
                @test rs2 == rs && rs2 !== rs
                @test issubset(rs2, rs)
                @test issubset(rs, rs2)
                delete!(rs2, 0)
                @test issubset(rs2, rs)
                @test !issubset(rs, rs2)

                rs = $TypeURSS{$Ti}((0:0, 2:4))
                @test filter(s->length(s)==3, rs) == $TypeURSS{$Ti}((2:4))
                @test (filter!(s->length(s)==1, rs); rs == $TypeURSS{$Ti}((0:0)))

            end
        end
    end
end


@testset "show" begin
    for Ti in (Int, UInt, UInt16, Float64)
        for TypeURSS in (UnitRangesSortedVector, UnitRangesSortedSet)
            @eval begin

                io = IOBuffer()
                rs = $TypeURSS{$Ti}((0:0, 2:4))
                show(io, rs)
                @test String(take!(io)) == "$TypeURSS{$Ti}(" *
                                           string($Ti(0)) * ":" * string($Ti(0)) * ", " *
                                           string($Ti(2)) * ":" * string($Ti(4)) * ")"

                io = IOBuffer()
                rs = $TypeURSS{$Ti}((0:0, 2:4))
                show(io, MIME("text/plain"), rs)
                @test String(take!(io)) == "$TypeURSS{$Ti}():\n" *
                                           "  " * string($Ti(0)) * ":" * string($Ti(0)) * "\n" *
                                           "  " * string($Ti(2)) * ":" * string($Ti(4))

            end
        end
    end
end


