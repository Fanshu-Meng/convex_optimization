# 定义 Momentum 优化算法函数
function momentum(f, ∇f, ∇²f, x0; α = 0.01, γ = 0.9, max_iter = 1000, tol = 1e-6)
    # 初始化参数
    x = copy(x0)
    # 初始化动量
    v = zeros(length(x))
    # 初始化迭代次数
    iter = 0
    # 初始化收敛曲线
    convergence_curve = Float64[]
    # 初始化迭代点列表
    iter_points = [copy(x)]

    # 开始迭代
    while iter < max_iter
        # 计算梯度
        grad = ∇f(x)
        # 更新动量
        v = γ * v + α * grad
        # 更新参数
        x = x - v

        # 记录函数值
        push!(convergence_curve, f(x))
        # 记录迭代点
        push!(iter_points, copy(x))

        # 判断是否收敛
        if norm(grad) < tol
            break
        end

        # 增加迭代次数
        iter += 1
    end

    return x, iter_points, convergence_curve
end