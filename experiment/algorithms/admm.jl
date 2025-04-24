function admm(f, ∇f, ∇²f, x0; ρ=1.0, max_iter=100, tol=1e-6)
    n = length(x0)
    # 初始化主变量和辅助变量
    x = copy(x0)
    z = copy(x0)
    # 初始化缩放后的对偶变量
    u = zeros(n)

    # 保存迭代点和函数值曲线
    iter_points = [copy(x)]
    convergence_curve = [f(x)]

    for k in 1:max_iter
        # x 更新：解线性方程 (∇²f(x) + ρ I) dx = - (∇f(x) + ρ (x - z + u))
        H = ∇²f(x)
        A = H + ρ * I(n)
        b = - (∇f(x) .+ ρ * (x .- z .+ u))
        dx = A \ b
        x = x .+ dx

        # z 更新：z = x + u
        z_old = copy(z)
        z = x .+ u

        # u 更新：u = u + (x - z)
        u = u .+ (x .- z)

        # 记录迭代信息
        push!(iter_points, copy(x))
        push!(convergence_curve, f(x))

        # 判断收敛性：原始残差和对偶残差
        r_norm = norm(x .- z)
        s_norm = norm(-ρ * (z .- z_old))
        if r_norm < tol && s_norm < tol
            # println("ADMM 在第 $k 次迭代收敛 (r_norm=$(round(r_norm, sigdigits=4)), s_norm=$(round(s_norm, sigdigits=4)))")
            break
        end
    end

    return x, iter_points, convergence_curve
end