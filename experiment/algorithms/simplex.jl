# Nelder-Mead 单纯形法实现（适用于无约束优化）

using Statistics

function simplex(f, ∇f, ∇²f, x0; α=1.0, γ=2.0, ρ=0.5, σ=0.5, max_iter=200, tol=1e-6)
    n = length(x0)
    # 初始化单纯形
    simplex = [x0]
    for i in 1:n
        ei = zeros(n)
        ei[i] = 0.05
        push!(simplex, x0 + ei)
    end

    iter_points = [x0]
    convergence_curve = [f(x0)]

    for iter in 1:max_iter
        # 按函数值排序
        sort!(simplex, by=f)
        best = simplex[1]
        worst = simplex[end]
        second_worst = simplex[end-1]

        # 计算质心（除最差点外）
        centroid = sum(simplex[1:end-1]) / n

        # 反射
        xr = centroid + α * (centroid - worst)
        fr = f(xr)

        if f(best) <= fr < f(second_worst)
            simplex[end] = xr
        elseif fr < f(best)
            # 扩展
            xe = centroid + γ * (xr - centroid)
            if f(xe) < fr
                simplex[end] = xe
            else
                simplex[end] = xr
            end
        else
            # 收缩
            xc = centroid + ρ * (worst - centroid)
            if f(xc) < f(worst)
                simplex[end] = xc
            else
                # 缩小单纯形
                for i in 2:length(simplex)
                    simplex[i] = best + σ * (simplex[i] - best)
                end
            end
        end

        # 记录路径
        push!(iter_points, simplex[1])
        push!(convergence_curve, f(simplex[1]))

        # 收敛判断
        if std([f(x) for x in simplex]) < tol
            break
        end
    end

    return simplex[1], iter_points, convergence_curve
end
