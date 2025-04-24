function subgradient_descent(f, ∇f, ∇²f, x0; max_iter=1000, tol=1e-6, α=0.01)
    # 初始化当前点为初始点
    x = copy(x0)
    # 存储迭代点的数组
    iter_points = [copy(x)]
    # 存储每次迭代的函数值
    convergence_curve = [f(x)]
    # 开始迭代
    for i in 1:max_iter
        # 计算当前点的子梯度
        subgrad = ∇f(x)
        # 更新当前点
        x = x - α * subgrad
        # 将更新后的点添加到迭代点数组中
        push!(iter_points, copy(x))
        # 计算更新后的函数值
        f_val = f(x)
        # 将更新后的函数值添加到收敛曲线数组中
        push!(convergence_curve, f_val)
        # 判断是否满足收敛条件
        if norm(subgrad) < tol
            break
        end
    end
    # 返回最优解、迭代点数组和收敛曲线数组
    return x, iter_points, convergence_curve
end