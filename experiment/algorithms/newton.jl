using LinearAlgebra

# 牛顿法
function newton(f, ∇f, ∇²f, x0; max_iter = 1000, tol = 1e-6)
    x = copy(x0)
    iter_points = [copy(x)]
    convergence_curve = [f(x)]
    for i in 1:max_iter
        g = ∇f(x)
        H = ∇²f(x)
        p = -H \ g
        x = x + p
        push!(iter_points, copy(x))
        push!(convergence_curve, f(x))
        if norm(g) < tol
            break
        end
    end
    return x, iter_points, convergence_curve
end    