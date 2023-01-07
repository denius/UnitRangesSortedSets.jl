
using BenchmarkTools
using DataFrames
using DataStructures
using PrettyTables
using Random
using UnitRangesSortedSets
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

function testfun_create_ext(T::Type, idx)
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

function testfun_in_ext(rs, idx)
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


function output_results(results)
    println("\nSizes:")
    PrettyTables.pretty_table(results[:sz         ], tf=PrettyTables.tf_markdown)
    println("\nIterate ranges consecutively, ms:")
    PrettyTables.pretty_table(results[:iter_ranges], tf=PrettyTables.tf_markdown)
    println("\nIterate elements consecutively, ms:")
    PrettyTables.pretty_table(results[:iter_con   ], tf=PrettyTables.tf_markdown)
    println("\nFill element-wise consecutively, ms:")
    PrettyTables.pretty_table(results[:fill_con   ], tf=PrettyTables.tf_markdown)
    println("\nAccess element-wise consecutively (time per element, ns):")
    PrettyTables.pretty_table(results[:accss_con  ], tf=PrettyTables.tf_markdown)
    println("\nFill element-wise in random order of elements, ms:")
    PrettyTables.pretty_table(results[:fill_rnd   ], tf=PrettyTables.tf_markdown)
    println("\nAccess element-wise randomly (time per element, ns):")
    PrettyTables.pretty_table(results[:accss_rnd  ], tf=PrettyTables.tf_markdown)
    nothing
end


function run_bench(name, N, indices, results, verbose = false)

    indices = sort(indices)

    Random.seed!(1234)
    rnd_indices = shuffle(indices)

    sizes = []
    push!(sizes, name)
    push!(sizes, length(indices))
    push!(sizes, length(testfun_create_ext(UnitRangesSortedSet{Int}, indices)))

    println("\nSizes:")

    for T in (BitSet, Set{Int}, SortedSet{Int}, UnitRangesSortedSet{Int}, VecUnitRangesSortedSet{Int})
        push!(sizes, Base.summarysize( testfun_create_ext(T, indices) ) )
    end
    push!(results[:sz], sizes)
    verbose && PrettyTables.pretty_table(results[:sz], tf=PrettyTables.tf_markdown)

    times_ir,
    times_ic,
    times_fc,
    times_ac,
    times_fr,
    times_ar = ntuple(_->Any[name, length(indices), length(testfun_create_ext(UnitRangesSortedSet{Int}, indices))], 6)

    for T in (BitSet, Set{Int}, SortedSet{Int}, UnitRangesSortedSet{Int}, VecUnitRangesSortedSet{Int})
        @sync begin
            # Iterate ranges consecutively, ms
            Threads.@spawn push!(times_ir, round(sigdigits=3, mean( @benchmark testfun_iter_ranges(r, $N) setup = (r = testfun_create_ext($T, $indices) ) ; ).time / 1e6))
            # Iterate elements consecutively, ms
            Threads.@spawn push!(times_ic, round(sigdigits=3, mean( @benchmark testfun_iter_cons(r, $N) setup = (r = testfun_create_ext($T, $indices) ) ; ).time / 1e6))
            # Fill element-wise consecutively, ms
            Threads.@spawn push!(times_fc, round(sigdigits=3, mean( @benchmark testfun_create_ext($T, $indices)                        ; ).time / 1e6))
            # Access element-wise consecutively (time per element, ns
            Threads.@spawn push!(times_ac, round(sigdigits=3, mean( @benchmark testfun_in_ext(r, $indices) setup = (r = testfun_create_ext($T, $indices) ) ; ).time / N))
            # Fill element-wise in random order of elements, ms
            Threads.@spawn push!(times_fr, round(sigdigits=3, mean( @benchmark testfun_create_ext($T, $rnd_indices)                        ; ).time / 1e6))
            # Access element-wise randomly (time per element, ns)
            Threads.@spawn push!(times_ar, round(sigdigits=3, mean( @benchmark testfun_in_ext(r, $rnd_indices) setup = (r = testfun_create_ext($T, $indices) ) ; ).time / N))
        end
    end
    push!(results[:iter_ranges], times_ir)
    push!(results[:iter_con], times_ic)
    push!(results[:fill_con], times_fc)
    push!(results[:accss_con], times_ac)
    push!(results[:fill_rnd], times_fr)
    push!(results[:accss_rnd], times_ar)

    verbose && output_results(results)

end

function run_benches(N = 1_000_000; verbose = false)

    # Size, bytes => Int,
    sz = DataFrame("ρ, %" => Float64[],
                   "nelems" => Int[],
                   "nranges" => Int[],
                   "BitSet" => Int[],
                   "Set" => Int[],
                   "SortedSet" => Int[],
                   "URSSet" => Int[],
                   "VURSSet" => Int[]
                  )


    iter_ranges,
    iter_con,
    fill_con,
    accss_con,
    fill_rnd,
    accss_rnd = ntuple(_->DataFrame("ρ, %" => Float64[],
                                    "nelems" => Int[],
                                    "nranges" => Int[],
                                    "BitSet" => Float64[],
                                    "Set" => Float64[],
                                    "SortedSet" => Float64[],
                                    "URSSet" => Float64[],
                                    "VURSSet" => Float64[]
                                   )
                       , 6)

    results = Dict((
                   :sz          => sz,
                   :iter_ranges => iter_ranges,
                   :iter_con    => iter_con,
                   :fill_con    => fill_con,
                   :accss_con   => accss_con,
                   :fill_rnd    => fill_rnd,
                   :accss_rnd   => accss_rnd
                  ))



    for density in (0.001, 0.01, 0.1, 1.0, 10.0, 50.0, 90., 99.0, 99.9, 99.99, 99.999)

        fraction = density / 100.0

        println("\n\ndensity = ", density, "%")

        nchunks = 1000
        chunklen = max(1, floor(Int, N / nchunks))
        Random.seed!(1234)
        indices = Int[]
        for i = 0:nchunks-1
            len = floor(Int, chunklen*fraction + randn() * chunklen * min(0.1, (1.0-fraction), fraction))
            len = max(1, min(chunklen-2, len))
            for j = 1:len
                push!(indices, i*chunklen + j)
            end
        end

        indices = unique(sort(indices))

        indices[1] != 1 && prepend!(indices, 1)
        indices[end] != N && append!(indices, N)

        run_bench(density, N, indices, results, verbose)

    end

    output_results(results)

    for density in (0.001, 0.01, 0.1, 1.0, 10.0, 50.0, 90., 99.0, 99.9, 99.99, 99.999)

        println("\n\ndensity = ", density, "%")

        Random.seed!(1234)
        indices = randsubseq(1:N, density / 100.0)
        indices[1] != 1 && prepend!(indices, 1)
        indices[end] != N && append!(indices, N)

        run_bench(density, N, indices, results, verbose)

    end

    output_results(results)

end


