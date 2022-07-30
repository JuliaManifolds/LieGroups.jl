function skewsymmetric(x::AbstractVector)
    N = length(x)
    if N == 2
        E1 = [0 -1;
              1  0]
        return x[1]*E1
    elseif N == 3
        E1 = [0 0  0;
              0 0 -1;
              0 1  0]
        E2 = [ 0 0 1;
               0 0 0;
              -1 0 0]
        E3 = [0 -1 0;
              1  0 0;
              0  0 0]
        return x[1]*E1 + x[2]*E2 + x[3]*E3
    else
        throw(ArgumentError("not support."))
    end
end

is_dof(::Type{T}, d::Int) where {T} = dof(T) == d
is_dim(::Type{T}, d::Int) where {T} = dim(T) == d
isskewsymmetric(A::AbstractMatrix) = A' == -A

check_dof(::Type{T}, d::Int) where {T} = @assert is_dof(T, d)
check_dim(::Type{T}, d::Int) where {T} = @assert is_dim(T, d)
check_skewsymmetric(X::AbstractMatrix)  = @assert isskewsymmetric(X)

euler_axis(R::AbstractArray) = real.(eigvecs(R)[:, size(R, 1)])

rotation_angle(R::AbstractArray) = acos(0.5 * (tr(R) - 1.))
