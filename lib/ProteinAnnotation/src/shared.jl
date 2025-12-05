using OneHotArrays

struct Categorizer{T}
    cats::T
    max::Int
    default::Int
end

Categorizer(c,m) = Categorizer(c,m,0)
Categorizer(a::AbstractVector{<:Real}) = Categorizer(v -> bin(v, a), length(a)+1)

bin(v, a::AbstractVector) = searchsortedfirst(vcat(a, [Inf]), v)

categorize(v::Union{Nothing,Missing}, d::Categorizer) = d.default
categorize(v::Union{Number,String,Char,Symbol}, d::Categorizer{<:Dict}) = get(d.cats, v, d.default)
categorize(v::Union{Number,String,Char,Symbol}, d::Categorizer{<:Function}) = d.cats(v)
categorize(v::AbstractArray, d::Categorizer) = categorize.(v, (d,))

onehotcats(v, d::Categorizer) = onehotbatch(categorize(v, d), 0:d.max) .* true
