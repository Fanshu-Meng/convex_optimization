# 三次正则化牛顿法
function cubic_regularized_newton(f, ∇f, ∇²f, x0; max_iter = 100, tol = 1e-6, γ = 1.0)
    x = copy(x0)
    iter_points = [copy(x)]
    convergence_curve = [f(x)]

    for i in 1:max_iter
        g = ∇f(x)
        H = ∇²f(x)

        # 求解子问题以获取搜索方向
        function cubic_subproblem(H, g, γ)
            # 这里可以使用更高效的方法求解子问题
            n = length(g)
            p = zeros(n)
            λ = 0.0
            while true
                try
                    p = - (H + λ * I) \ g
                    break
                catch
                    λ += 1e-3
                end
            end
            return p
        end

        p = cubic_subproblem(H, g, γ)

        # 更新点
        x = x + p

        push!(iter_points, copy(x))
        push!(convergence_curve, f(x))

        # 检查收敛条件
        if norm(p) < tol
            break
        end
    end

    return x, iter_points, convergence_curve
end