# Krylov 子空间法用于无约束优化
function krylov(f, ∇f, ∇²f, x0; max_iter=100, tol=1e-6)
    x = x0
    fx = f(x)
    grad = ∇f(x)
    H = ∇²f(x)

    iter_points = [x]
    convergence_curve = [fx]

    for iter in 1:max_iter
        # 判断梯度是否足够小，满足停止准则
        if norm(grad) < tol
            break
        end

        # 使用共轭梯度法求解 H * p = -grad，即在Krylov子空间内找搜索方向
        p = conjugate_gradient_solver(H, -grad)

        # 使用单位步长（可改进为线搜索）
        x = x + p

        # 更新函数值和梯度
        fx = f(x)
        grad = ∇f(x)
        H = ∇²f(x)

        # 记录迭代路径和收敛曲线
        push!(iter_points, x)
        push!(convergence_curve, fx)
    end

    return x, iter_points, convergence_curve
end

# 共轭梯度法求解线性系统 A * x = b
function conjugate_gradient_solver(A, b; tol=1e-8, max_iter=100)
    x = zeros(length(b))
    r = b - A * x
    p = copy(r)
    rs_old = dot(r, r)

    for i in 1:max_iter
        Ap = A * p
        alpha = rs_old / dot(p, Ap)
        x += alpha * p
        r -= alpha * Ap

        if norm(r) < tol
            break
        end

        rs_new = dot(r, r)
        beta = rs_new / rs_old
        p = r + beta * p
        rs_old = rs_new
    end

    return x
end
