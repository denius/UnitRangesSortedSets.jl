
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

#=
Sizes:
| Density, % | # elements | # ranges | BitSet |      Set | SortedSet |   URSSet | VURSSet |
|    Float64 |      Int64 |    Int64 |  Int64 |    Int64 |     Int64 |    Int64 |   Int64 |
|------------|------------|----------|--------|----------|-----------|----------|---------|
|      0.001 |         13 |       13 | 125064 |      768 |      1200 |     1336 |     312 |
|       0.01 |         90 |       90 | 125064 |     2496 |      6136 |     6888 |    1544 |
|        0.1 |        985 |      985 | 125064 |    37056 |     63336 |    71248 |   15864 |
|        1.0 |      10126 |    10014 | 125064 |   147648 |    649600 |   722608 |  160328 |
|       10.0 |      99784 |    89966 | 125064 |  2359488 |   6398824 |  6488904 | 1439560 |
|       50.0 |     499795 |   249673 | 125064 |  9437376 |  32049480 | 18007680 | 3994872 |
|       90.0 |     900008 |    89882 | 125064 | 18874560 |  57712944 |  6482848 | 1438216 |
|       99.0 |     989970 |     9934 | 125064 | 18874560 |  63481904 |   716792 |  159048 |
|       99.9 |     998975 |     1025 | 125064 | 18874560 |  64059304 |    74376 |   16504 |
|      99.99 |     999898 |      103 | 125064 | 18874560 |  64118488 |     7824 |    1752 |
|     99.999 |     999985 |       16 | 125064 | 18874560 |  64124112 |     1648 |     360 |

Iterate ranges consecutively, ms:
| Density, % |     BitSet |     Set |   SortedSet |      URSSet |     VURSSet |
|    Float64 |    Float64 | Float64 |     Float64 |     Float64 |     Float64 |
|------------|------------|---------|-------------|-------------|-------------|
|      0.001 | 0.00470375 | 32.9019 | 0.000237866 | 0.000228618 |  2.10459e-5 |
|       0.01 |  0.0055848 | 21.8019 |  0.00144485 |   0.0014535 |  3.20422e-5 |
|        0.1 |  0.0149049 |  34.726 |   0.0166525 |   0.0172187 |  3.91753e-5 |
|        1.0 |  0.0715772 | 53.2173 |    0.170203 |    0.176256 | 0.000122602 |
|       10.0 |   0.183567 | 44.1411 |     2.54265 |     2.29348 | 0.000758026 |
|       50.0 |   0.572195 | 48.4386 |     11.7606 |     6.07555 |  0.00207106 |
|       90.0 |   0.947723 | 29.2782 |     22.4078 |     2.28402 |  0.00075313 |
|       99.0 |     1.0006 | 24.0011 |     24.4706 |    0.183785 | 0.000113183 |
|       99.9 |   0.994607 | 24.9677 |     26.4903 |   0.0228217 |  4.34554e-5 |
|      99.99 |   0.974751 | 26.4434 |     25.7705 |  0.00167289 |  2.39065e-5 |
|     99.999 |    1.00035 | 26.0407 |     25.2084 | 0.000272213 |  2.55414e-5 |

Iterate elements consecutively, ms:
| Density, % |    BitSet |     Set |   SortedSet |      URSSet |     VURSSet |
|    Float64 |   Float64 | Float64 |     Float64 |     Float64 |     Float64 |
|------------|-----------|---------|-------------|-------------|-------------|
|      0.001 | 0.0123706 | 9.14295 | 0.000209997 | 0.000239995 |  3.80409e-5 |
|       0.01 | 0.0128422 | 16.1447 |  0.00129227 | 0.000756997 | 0.000122229 |
|        0.1 | 0.0177811 | 14.3266 |   0.0152082 |   0.0182257 |  0.00115038 |
|        1.0 | 0.0902567 | 24.2161 |    0.161449 |     0.18974 |   0.0114765 |
|       10.0 |  0.144444 | 22.1168 |     2.44871 |     2.40601 |    0.112726 |
|       50.0 |  0.259572 | 31.5915 |     11.8219 |     6.27201 |    0.341234 |
|       90.0 |   0.36913 | 25.1313 |     21.0716 |     2.43445 |    0.116856 |
|       99.0 |  0.347011 | 24.2367 |     23.5105 |    0.200395 |   0.0128309 |
|       99.9 |  0.311916 | 25.3368 |     24.0629 |   0.0247403 |  0.00125932 |
|      99.99 |  0.287967 | 24.9056 |     24.1566 |  0.00183615 | 0.000148431 |
|     99.999 |  0.285591 | 23.7365 |     24.0033 |  0.00029194 |  4.13598e-5 |

Fill element-wise consecutively, ms:
| Density, % |    BitSet |        Set |  SortedSet |     URSSet |    VURSSet |
|    Float64 |   Float64 |    Float64 |    Float64 |    Float64 |    Float64 |
|------------|-----------|------------|------------|------------|------------|
|      0.001 | 0.0504952 |  0.0010947 | 0.00215711 | 0.00263551 | 0.00112078 |
|       0.01 | 0.0283609 | 0.00332908 | 0.00743078 |  0.0120598 | 0.00543543 |
|        0.1 | 0.0413336 |  0.0268798 |   0.106045 |   0.194139 |  0.0489959 |
|        1.0 |  0.259884 |   0.280313 |    1.40713 |    2.26005 |    0.50943 |
|       10.0 |  0.743347 |    3.17031 |    13.9822 |    22.9465 |    4.87803 |
|       50.0 |   1.66972 |    22.2561 |    71.8956 |    84.6329 |    22.1698 |
|       90.0 |    2.2679 |    46.1752 |    144.826 |    84.7104 |    25.1256 |
|       99.0 |   3.75522 |    50.5247 |     159.61 |    69.9275 |    21.9703 |
|       99.9 |    3.1204 |    52.8201 |    180.795 |    61.6057 |    21.3702 |
|      99.99 |   2.71726 |    49.8933 |    154.492 |    53.1448 |    18.3463 |
|     99.999 |    2.7054 |    47.3519 |    159.132 |    48.3507 |    17.1256 |

Access element-wise consecutively (time per element, ns):
| Density, % |      BitSet |         Set |   SortedSet |      URSSet |     VURSSet |
|    Float64 |     Float64 |     Float64 |     Float64 |     Float64 |     Float64 |
|------------|-------------|-------------|-------------|-------------|-------------|
|      0.001 |  3.92055e-5 |  0.00014986 |  7.41426e-5 | 0.000141601 |  8.99724e-5 |
|       0.01 |  7.29371e-5 | 0.000935166 | 0.000760177 |  0.00141584 | 0.000633098 |
|        0.1 | 0.000628493 |   0.0109913 |   0.0155518 |   0.0234083 |   0.0290949 |
|        1.0 |  0.00494196 |    0.166526 |     0.25676 |    0.332913 |     0.36437 |
|       10.0 |   0.0480806 |     1.64649 |      3.8328 |     4.52884 |     3.71688 |
|       50.0 |    0.246509 |     10.0217 |     20.9038 |     17.6052 |     12.4291 |
|       90.0 |    0.445322 |     21.0237 |      44.603 |     8.27193 |     5.23677 |
|       99.0 |    0.479589 |     23.4211 |     49.6495 |     3.36669 |     1.67657 |
|       99.9 |    0.510899 |     25.0122 |     50.6898 |     2.97167 |     1.39197 |
|      99.99 |     0.48278 |     24.3924 |     49.7275 |     2.94061 |     1.38486 |
|     99.999 |    0.507975 |     22.3382 |     51.1434 |     2.90628 |     1.29976 |

Fill element-wise in random order of elements, ms:
| Density, % |    BitSet |         Set |  SortedSet |     URSSet |    VURSSet |
|    Float64 |   Float64 |     Float64 |    Float64 |    Float64 |    Float64 |
|------------|-----------|-------------|------------|------------|------------|
|      0.001 | 0.0297153 | 0.000550401 | 0.00194611 | 0.00266436 | 0.00142334 |
|       0.01 | 0.0123559 |  0.00301805 | 0.00690287 |  0.0118373 | 0.00821339 |
|        0.1 | 0.0518161 |   0.0306109 |   0.103354 |   0.215148 |   0.143558 |
|        1.0 | 0.0394055 |    0.281552 |    1.47494 |     2.6985 |    4.07263 |
|       10.0 |   0.24699 |     3.07479 |    23.9683 |    38.5171 |    271.336 |
|       50.0 |   1.04152 |     20.1027 |    160.924 |    207.269 |    6053.56 |
|       90.0 |   1.79929 |     43.9382 |    361.069 |    367.019 |    12259.5 |
|       99.0 |   2.35568 |     48.4441 |    446.941 |    417.042 |    12596.6 |
|       99.9 |   1.91882 |     47.8206 |    390.915 |    391.937 |    12324.1 |
|      99.99 |   1.84786 |     51.2692 |    417.785 |    405.034 |    12456.1 |
|     99.999 |   2.51904 |     52.6851 |    413.885 |      410.4 |    12517.1 |

Access element-wise randomly (time per element, ns):
| Density, % |      BitSet |         Set |   SortedSet |      URSSet |     VURSSet |
|    Float64 |     Float64 |     Float64 |     Float64 |     Float64 |     Float64 |
|------------|-------------|-------------|-------------|-------------|-------------|
|      0.001 |   3.8967e-5 | 0.000150629 |  7.36069e-5 | 0.000143181 |  9.05463e-5 |
|       0.01 |  7.24929e-5 | 0.000950753 | 0.000623426 |  0.00117391 | 0.000729448 |
|        0.1 | 0.000556236 |   0.0112985 |   0.0362643 |   0.0579239 |   0.0489833 |
|        1.0 |  0.00542677 |    0.138121 |    0.713657 |    0.906538 |    0.674067 |
|       10.0 |   0.0493276 |     1.60159 |     12.7708 |     15.9306 |     8.85079 |
|       50.0 |    0.287777 |     7.71945 |     91.1639 |     99.4458 |     53.7053 |
|       90.0 |     0.51347 |     20.4769 |      217.11 |     149.032 |     75.2164 |
|       99.0 |    0.557918 |     24.2201 |     275.467 |     99.3213 |     65.5492 |
|       99.9 |    0.558514 |     23.4907 |     263.965 |     56.8168 |     48.5415 |
|      99.99 |    0.543366 |     23.6841 |     262.209 |     32.5893 |     31.0012 |
|     99.999 |    0.572952 |     23.4482 |     277.096 |     19.8198 |     17.8942 |
=#

