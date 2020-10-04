struct InternalArray{T,N,P} <: AbstractArray{T,N}
	parent::P
	extracted::Int
end

#TODO check the dimensions and everything is as expected
# That is: parent[i][extracted] is meaningful
function InternalArray(parent::AbstractArray, extracted::Int)
	InternalArray{eltype(parent), ndims(eltype(parent)), typeof(parent)}(parent, extracted)
end

Base.size(v::InternalArray) = size(v.parent)
Base.getindex(v::InternalArray, i) = v.parent[i][v.extracted]

