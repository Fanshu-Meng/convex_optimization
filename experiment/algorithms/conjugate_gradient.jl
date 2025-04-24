using LinearAlgebra

# 共轭梯度法
function conjugate_gradient(f, ∇f, ∇²f, x0; max_iter = 1000, tol = 1e-6)
    x = copy(x0)
    iter_points = [copy(x)]
    convergence_curve = [f(x)]
    g = ∇f(x)
    d = -g
    for i in 1:max_iter
        α = line_search(f, ∇f, x, d)
        x = x + α * d
        push!(iter_points, copy(x))
        push!(convergence_curve, f(x))
        g_new = ∇f(x)
        β = dot(g_new, g_new - g) / dot(g, g)
        d = -g_new + β * d
        g = g_new
        if norm(g) < tol
            break
        end
    end
    return x, iter_points, convergence_curve
end

# 简单的线搜索函数
function line_search(f, ∇f, x, d)
    α = 1.0
    ρ = 0.5
    c = 0.1
    while f(x + α * d) > f(x) + c * α * dot(∇f(x), d)
        α = ρ * α
    end
    return α
end    