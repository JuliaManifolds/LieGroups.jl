isskewsymmetric(A::AbstractMatrix) = A' == -A

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

euler_axis(R::AbstractArray) = real.(eigvecs(R)[:, size(R, 1)])

rotation_angle(R::AbstractArray) = acos(0.5 * (tr(R) - 1.))
