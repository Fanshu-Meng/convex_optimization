using LinearAlgebra

# 阻尼牛顿法
function damped_newton(f, ∇f, ∇²f, x0; max_iter = 1000, tol = 1e-6)
    x = copy(x0)
    iter_points = [copy(x)]
    convergence_curve = [f(x)]
    for i in 1:max_iter
        g = ∇f(x)
        H = ∇²f(x)
        p = -H \ g
        α = line_search(f, ∇f, x, p)
        x = x + α * p
        push!(iter_points, copy(x))
        push!(convergence_curve, f(x))
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