using LinearAlgebra

# 共轭方向法实现
function conjugate_direction(f, ∇f, ∇²f, x0; max_iter=1000, tol=1e-6)
    # 初始化当前点
    x = copy(x0)
    # 初始化迭代次数
    iter = 0
    # 初始化迭代点列表
    iter_points = [copy(x)]
    # 初始化收敛曲线列表
    convergence_curve = [f(x)]
    # 获取变量的维度
    n = length(x)
    # 初始化共轭方向集
    directions = [zeros(n) for _ in 1:n]
    # 初始搜索方向为负梯度方向
    d = -∇f(x)

    while iter < max_iter
        # 计算步长，这里简单使用固定步长 0.1，可根据需要修改
        α = line_search(f, ∇f, x, d)
        # 更新当前点
        x = x + α * d
        # 记录当前迭代点
        push!(iter_points, copy(x))
        # 记录当前函数值
        push!(convergence_curve, f(x))

        # 计算新的负梯度方向
        new_d = -∇f(x)
        # 计算 β 值，这里使用 Fletcher-Reeves 公式
        β = dot(new_d, new_d) / dot(d, d)
        # 更新搜索方向
        d = new_d + β * d

        # 判断是否满足收敛条件
        if norm(new_d) < tol
            break
        end

        # 更新迭代次数
        iter += 1
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