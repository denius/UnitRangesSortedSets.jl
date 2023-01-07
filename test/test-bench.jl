
using BenchmarkTools
using DataFrames
using DataStructures
using PrettyTables
using Random
using .UnitRangesSortedSets
using Unitful

#
#  Testing functions
#

function testfun_create(T::Type, N = 1_000_000, density = 0.9)
    rs = T()
    Random.seed!(1234)
    randseq = randsubseq(1:N, density)
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

function testfun_create_cons(T::Type, N = 1_000_000, density = 0.9)
    rs = T()
    Random.seed!(1234)
    randseq = randsubseq(1:N, density)
    for i in randseq
        push!(rs, i)
    end
    for i in randseq
        in(i, rs) || println("Not coincide on index $i")
    end
    sum(length(r) for r in rs) == length(randseq) || println("Lost some indices")
    rs
end

function testfun_create_outer(T::Type, idx)
    rs = T()
    for i in idx
        push!(rs, i)
    end
    rs
end

function testfun_create_dense(T::Type, N = 1_000_000, nchunks = 800, density = 0.95)
    rs = T()
    chunklen = max(1, floor(Int, N / nchunks))
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
    num = 0
    for i in indices
        in(i, rs) && (num += 1)
    end
    num
end

function testfun_in_outer(rs, idx)
    num = 0
    for i in idx
        in(i, rs) && (num += 1)
    end
    num
end

function testfun_in_rand(rs)
    len = sum(length(r) for r in rs)
    Random.seed!(1234)
    num = 0
    for j = 1:len
        i = rand(1:len)
        in(i, rs) && (num += 1)
    end
    num
end

function testfun_in_cons(rs)
    start = first(first(rs))
    stop = last(last(rs))
    num = 0
    for i in start:stop
        in(i, rs) && (num += 1)
    end
    num
end

function testfun_iter_cons(rs::AbstractUnitRangesSortedSet, N = length(rs))
    num = 0
    for r in rs, i in r
        num += 1
    end
    num
end

function testfun_iter_cons(rs::Union{BitSet,SortedSet}, N = 0)
    num = 0
    for i in rs
        num += 1
    end
    num
end

function testfun_iter_cons(rs::Set, N)
    num = 0
    inrange = false
    for i = 1:N
        in(i, rs) && (num += 1)
    end
    num
end

function testfun_iter_ranges(rs::AbstractUnitRangesSortedSet, N = length(rs))
    num = 0
    for r in rs
        num += 1
    end
    num
end

function testfun_iter_ranges(rs::Union{BitSet,SortedSet}, N = 0)
    num = 0
    prev = first(rs)
    for i in rs
        if i != prev + 1
            num += 1
        end
        prev = i
    end
    num
end
function testfun_iter_ranges(rs::Set, N)
    num = 0
    inrange = false
    for i = 1:N
        if in(i, rs) && inrange
            continue
        elseif in(i, rs) && !inrange
            inrange = true
            num += 1
        elseif !in(i, rs) && inrange
            inrange = false
        else # !in(i, rs) && !inrange
            continue
        end
    end
    num
end
#function testfun_iter_ranges(rs::Union{BitSet,Set}, N)
#    num = 1
#    prev = 0
#    for i = 1:N
#        if i in rs
#            prev = i
#            break
#        end
#    end
#    for i = (prev+1):N
#        if i in rs && i != prev + 1
#            num += 1
#        end
#        prev = i
#    end
#    num
#end

function run_bench(N = 1_000_000)

    # Size, bytes => Int,
    sz = DataFrame("Density, %" => Float64[],
                   "# elements" => Int[],
                   "# ranges" => Int[],
                   "BitSet" => Int[],
                   "Set" => Int[],
                   "SortedSet" => Int[],
                   "URSSet" => Int[],
                   "VURSSet" => Int[]
                  )

    fill_con = DataFrame("Density, %" => Float64[],
                   "BitSet" => Float64[],
                   "Set" => Float64[],
                   "SortedSet" => Float64[],
                   "URSSet" => Float64[],
                   "VURSSet" => Float64[]
                  )

    fill_rnd = DataFrame("Density, %" => Float64[],
                   "BitSet" => Float64[],
                   "Set" => Float64[],
                   "SortedSet" => Float64[],
                   "URSSet" => Float64[],
                   "VURSSet" => Float64[]
                  )

    accss_con = DataFrame("Density, %" => Float64[],
                   "BitSet" => Float64[],
                   "Set" => Float64[],
                   "SortedSet" => Float64[],
                   "URSSet" => Float64[],
                   "VURSSet" => Float64[]
                  )

    accss_rnd = DataFrame("Density, %" => Float64[],
                   "BitSet" => Float64[],
                   "Set" => Float64[],
                   "SortedSet" => Float64[],
                   "URSSet" => Float64[],
                   "VURSSet" => Float64[]
                  )

    iter_ranges = DataFrame("Density, %" => Float64[],
                   "BitSet" => Float64[],
                   "Set" => Float64[],
                   "SortedSet" => Float64[],
                   "URSSet" => Float64[],
                   "VURSSet" => Float64[]
                  )

    iter_con = DataFrame("Density, %" => Float64[],
                   "BitSet" => Float64[],
                   "Set" => Float64[],
                   "SortedSet" => Float64[],
                   "URSSet" => Float64[],
                   "VURSSet" => Float64[]
                  )



    for density in (0.001, 0.01, 0.1, 1.0, 10.0, 50.0, 90., 99.0, 99.9, 99.99, 99.999)

        println("\n\ndensity = ", density, "%")

        Random.seed!(1234)
        indices = randsubseq(1:N, density / 100.0)
        indices[1] != 1 && prepend!(indices, 1)
        indices[end] != N && append!(indices, N)

        result = []
        push!(result, density)
        push!(result, length(indices))
        push!(result, length(testfun_create_outer(UnitRangesSortedSet{Int}, indices)))

        println("\nSizes:")

        for T in (BitSet, Set{Int}, SortedSet{Int}, UnitRangesSortedSet{Int}, VecUnitRangesSortedSet{Int})
            push!(result, Base.summarysize( testfun_create_outer(T, indices) ) )
        end
        push!(sz, result)
        PrettyTables.pretty_table(sz, tf=PrettyTables.tf_markdown)

        println("\nIterate ranges consecutively, ms:")
        times = Float64[density]
        for T in (BitSet, Set{Int}, SortedSet{Int}, UnitRangesSortedSet{Int}, VecUnitRangesSortedSet{Int})
            push!(times, mean( @benchmark testfun_iter_ranges(r, $N) setup = (r = testfun_create_outer($T, $indices) ) ; ).time / 1e6)
        end
        push!(iter_ranges, times)
        PrettyTables.pretty_table(iter_ranges, tf=PrettyTables.tf_markdown)

        println("\nIterate elements consecutively, ms:")
        times = Float64[density]
        for T in (BitSet, Set{Int}, SortedSet{Int}, UnitRangesSortedSet{Int}, VecUnitRangesSortedSet{Int})
            push!(times, mean( @benchmark testfun_iter_cons(r, $N) setup = (r = testfun_create_outer($T, $indices) ) ; ).time / 1e6)
        end
        push!(iter_con, times)
        PrettyTables.pretty_table(iter_con, tf=PrettyTables.tf_markdown)

        println("\nFill element-wise consecutively, ms:")
        times = Float64[density]
        for T in (BitSet, Set{Int}, SortedSet{Int}, UnitRangesSortedSet{Int}, VecUnitRangesSortedSet{Int})
            push!(times, mean( @benchmark testfun_create_outer($T, $indices)                        ; ).time / 1e6)
        end
        push!(fill_con, times)
        PrettyTables.pretty_table(fill_con, tf=PrettyTables.tf_markdown)

        println("\nAccess element-wise consecutively (time per element, ns):")
        times = Float64[density]
        for T in (BitSet, Set{Int}, SortedSet{Int}, UnitRangesSortedSet{Int}, VecUnitRangesSortedSet{Int})
            push!(times, mean( @benchmark  testfun_in_outer(r, $indices) setup = (r = testfun_create_outer($T, $indices) ) ; ).time / N)
        end
        push!(accss_con, times)
        PrettyTables.pretty_table(accss_con, tf=PrettyTables.tf_markdown)




        Random.seed!(1234)
        indices = shuffle(indices)


        println("\nFill element-wise in random order of elements, ms:")
        times = Float64[density]
        for T in (BitSet, Set{Int}, SortedSet{Int}, UnitRangesSortedSet{Int}, VecUnitRangesSortedSet{Int})
            push!(times, mean( @benchmark testfun_create_outer($T, $indices) ; ).time / 1e6)
        end
        push!(fill_rnd, times)
        PrettyTables.pretty_table(fill_rnd, tf=PrettyTables.tf_markdown)

        println("\nAccess element-wise randomly (time per element, ns):")
        Random.seed!(1234)
        shindices = shuffle(indices)
        times = Float64[density]
        for T in (BitSet, Set{Int}, SortedSet{Int}, UnitRangesSortedSet{Int}, VecUnitRangesSortedSet{Int})
            push!(times, mean( @benchmark  testfun_in_outer(r, $indices) setup = (r = testfun_create_outer($T, $indices) ) ; ).time / N)
        end
        push!(accss_rnd, times)
        PrettyTables.pretty_table(accss_rnd, tf=PrettyTables.tf_markdown)

    end

end

