# Experiment adapter for arXiv:2401.15492v1, not a runtime dependency.
# Includes unmodified pinned author files; assembly below follows the paper's
# row/column conventions explicitly rather than silently editing the README.
using LinearAlgebra, SparseArrays, SHA
BLAS.set_num_threads(1)

source = only(ARGS)
for (file, expected) in [
    ("collapse.jl", "4be26d62090ecce1c7fc1ff81033791931ee9ebf81b65af686cda52d7fb550ca"),
    ("cgls.jl", "e5c5d6b44065ea3d395fd6798520dd8b29ea2c84764b98ca4381d11f70aa2c61"),
    ("simplicialComplex.jl", "4e17b4872210e4f24a7d6c59973f38a41b9608fc0fe9057db455b6e2cb938e60"),
    ("utils.jl", "c8250f73c8794459f5673e6960742992fadfe98f391258ac228b4635d4c8ba16")]
    bytes2hex(sha256(read(joinpath(source, "src", file)))) == expected || error("unpinned author source: " * file)
end
include(joinpath(source, "src", "simplicialComplex.jl"))
include(joinpath(source, "src", "collapse.jl"))
include(joinpath(source, "src", "cgls.jl"))
include(joinpath(source, "src", "utils.jl"))

function json(value)
    if value isa AbstractString
        return "\"" * replace(value, "\\" => "\\\\", "\"" => "\\\"", "\n" => "\\n") * "\""
    elseif value isa AbstractDict
        return "{" * join((json(string(k)) * ":" * json(v) for (k, v) in sort(collect(value); by=first)), ",") * "}"
    elseif value isa AbstractArray
        return "[" * join(json.(value), ",") * "]"
    elseif value isa AbstractFloat && !isfinite(value)
        return "null"
    else
        return string(value)
    end
end

tokens = split(read(stdin, String))
position = 1
function takeint()
    global position
    value = parse(Int, tokens[position]); position += 1
    return value
end
function takefloat()
    global position
    value = parse(Float64, tokens[position]); position += 1
    return value
end
n, m, t, peel = takeint(), takeint(), takeint(), takeint()
edges = [takeint() + 1 for _ in 1:2, _ in 1:m] |> permutedims
faces = zeros(Int, t, 3)
weights = zeros(t)
for j in 1:t
    faces[j, :] = [takeint() + 1 for _ in 1:3]
    weights[j] = takefloat()
end
rhs = [takefloat() for _ in 1:m]
position == length(tokens) + 1 || error("trailing input")

if m == 0
    source_error = try
        greedyCollapse(edges, faces)
        "none"
    catch error
        sprint(showerror, error)
    end
    println(json(Dict("julia" => string(VERSION), "selected" => Int[],
        "pairs" => [], "selected_rank" => 0, "factor_error" => 0.,
        "status" => "empty_source_unsupported", "source_error" => source_error)))
    exit()
end

initial_ok, initial_edges, initial_faces, _, _, _, _ = greedyCollapse(edges, faces)
adjacency = getEdges2Trig(edges, faces)
triangle_edges = getTrig2Edg(edges, faces, adjacency)
selected = peel == 1 ? Int.(initial_faces) : Int[]
core = setdiff(collect(1:t), Int.(initial_faces))
candidates = peel == 1 ? core : collect(1:t)
sort!(candidates; by=i -> (-weights[i], i))
for triangle in candidates
    trial = [selected; triangle]
    edge_sets = [Set(intersect(collect(adj), trial)) for adj in adjacency]
    triangle_sets = [triangle_edges[j] for j in trial]
    ok, _, _, _, _, _, _ = greedyCollapseShort(edge_sets, triangle_sets, trial, edges)
    if ok
        push!(selected, triangle)
    end
end
local_faces = sort(selected)
ok, sigma, tau, _, _, _, _ = greedyCollapse(edges, faces[local_faces, :])
ok || error("selected complex did not collapse")
row_order = [Int.(sigma); setdiff(collect(1:m), Int.(sigma))]
column_order = local_faces[Int.(tau)]
r = length(column_order)
b = Matrix(B2fromTrig(edges, faces))
c = b[row_order, column_order] * Diagonal(sqrt.(weights[column_order]))
u = b * Diagonal(weights) * b'
up = u[row_order, row_order]
selected_u = b[:, selected] * Diagonal(weights[selected]) * b[:, selected]'
factor_error = norm(c * c' - selected_u[row_order, row_order]) / max(1., norm(selected_u))
selected_rank = r == 0 ? 0 : rank(b[:, selected])
full_rank = t == 0 ? 0 : rank(b)

# Retain diagnostics for the README's original acceptance-order lookup and
# Perm(v) convention. Both differ from direct row indexing on a general order.
literal_columns = selected[Int.(tau)]
literal_columns = [literal_columns; setdiff(collect(1:t), literal_columns)]
literal_c = Perm(row_order) * b * Diagonal(sqrt.(weights)) * indicator(selected, t) * Perm(literal_columns)
literal_upper_nnz = count(!iszero, triu(literal_c, 1))
paper_upper_nnz = count(!iszero, triu(c, 1))

status = "rank_loss"
original_residual = 0.
source_iterations = 0
source_x = Float64[]
active_spectrum = Float64[]
if selected_rank == full_rank
    cp = r == 0 ? zeros(0, m) : pinv(c; rtol=1e-12)
    a = cp * up * cp'
    h = cp * rhs[row_order]
    x, source_iterations = cgls(a, reshape(h, :, 1); tol=1e-12, maxit=4000)
    source_x = vec(x)
    if all(isfinite, x)
        lifted = zeros(m)
        lifted[row_order] = vec(cp' * x)
        original_residual = norm(u * lifted - rhs) / max(1., norm(rhs))
        status = "converged"
    else
        status = norm(rhs) == 0 ? "zero_rhs_source_nonfinite" : "source_nonfinite"
    end
    active_spectrum = eigvals(Symmetric(a))
end
println(json(Dict(
    "julia" => string(VERSION), "stdlib_only" => true,
    "selected" => selected .- 1,
    "pairs" => [[Int(e) - 1, local_faces[Int(j)] - 1] for (e, j) in zip(sigma, tau)],
    "initial_pairs" => [[Int(e) - 1, Int(j) - 1] for (e, j) in zip(initial_edges, initial_faces)],
    "selected_rank" => selected_rank, "full_rank" => full_rank,
    "factor_error" => factor_error, "paper_upper_nnz" => paper_upper_nnz,
    "literal_readme_upper_nnz" => literal_upper_nnz,
    "active_spectrum" => active_spectrum, "status" => status,
    "source_iterations" => source_iterations, "source_x" => source_x,
    "original_residual" => original_residual)))
