# 定义拟牛顿法函数
function quasi_newton(f, ∇f, ∇²f, x0; max_iter = 1000, tol = 1e-6)
    # 初始化迭代次数
    iter = 0
    # 初始化当前点为初始点
    x = x0
    # 初始化海森矩阵的近似为单位矩阵
    H = Matrix(I, length(x0), length(x0))
    # 存储迭代点
    iter_points = [x]
    # 存储每次迭代的函数值，用于收敛曲线
    convergence_curve = [f(x)]

    # 开始迭代
    while iter < max_iter
        # 计算当前点的梯度
        grad = ∇f(x)
        # 检查梯度的范数是否小于容差，如果是则认为收敛
        if norm(grad) < tol
            break
        end
        # 计算搜索方向，使用海森矩阵的近似进行计算
        d = -H * grad
        # 这里可以添加线搜索方法来确定步长，简单起见，这里假设步长为 1
        α = 1
        # 更新当前点
        x_new = x + α * d
        # 计算梯度的变化
        s = x_new - x
        # 计算梯度的变化
        y = ∇f(x_new) - grad
        # 使用 DFP 公式更新海森矩阵的近似
        ρ = 1 / (y' * s)[1]
        H = (I - ρ * s * y') * H * (I - ρ * y * s') + ρ * s * s'
        # 更新当前点
        x = x_new
        # 记录迭代点
        push!(iter_points, x)
        # 记录当前迭代的函数值
        push!(convergence_curve, f(x))
        # 迭代次数加 1
        iter += 1
    end
    return x, iter_points, convergence_curve
end    