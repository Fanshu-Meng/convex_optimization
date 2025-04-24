using SpecialFunctions

function DOA(f, ∇f, ∇²f, x0; max_iter=2000, pop_size=300, verbose=false)
    D = length(x0)  # 变量维度
    lower = -5.0 * ones(D)  # 搜索空间下界
    upper = 5.0 * ones(D)   # 搜索空间上界

    # 初始化种群
    population = [lower .+ rand(D) .* (upper .- lower) for _ in 1:pop_size]
    fitness = [f(x) for x in population]
    best_idx = argmin(fitness)
    best = copy(population[best_idx])
    
    iter_points = []  # 记录迭代路径
    convergence_curve = []  # 记录收敛曲线

    for iter in 1:max_iter
        alpha = rand() * ((iter / max_iter)^2 - 2 * iter / max_iter + 1)
        
        # 升空阶段（Rise stage）
        if randn() < 1.5
            theta = rand(pop_size) .* (2π) .- π
            y = randn(pop_size)
            y = broadcast((y_val) -> y_val < 0 ? 1e-6 : y_val, y)  # 处理负值，避免log报错
            InY = exp.(-0.5 * log.(y).^2) ./ (y .* sqrt(2π))
            InY = broadcast((y_val) -> y_val < 0 ? 0 : y_val, InY)  # 处理负值
            PopDec = [
                population[i] .+ alpha .* (cos(theta[i])/exp(theta[i])) .* (sin(theta[i])/exp(theta[i])) .* InY[i] .* (rand(D) .* (upper .- lower) .+ lower .- population[i])
                for i in 1:pop_size
            ]
        else
            PopDec = [
                population[i] .* (1 .- rand() * (((iter^2 - 2*iter + 1)/(max_iter^2 - 2*max_iter + 1)) + 1))
                for i in 1:pop_size
            ]
        end

        # 下降阶段（Decline stage）
        beta = randn(pop_size, D)
        mean_dec = mean(reduce(hcat, PopDec), dims=2)[:]
        PopDec = [
            PopDec[i] .- alpha .* beta[i, :] .* (mean_dec .- alpha .* beta[i, :] .* PopDec[i])
            for i in 1:pop_size
        ]

        # 着陆阶段（Land stage）
        fitness = [f(x) for x in PopDec]
        best_idx = argmin(fitness)
        best = copy(PopDec[best_idx])
        elite = [best for _ in 1:pop_size]
        C = (gamma(2.5) * sin(π * 0.75) / (gamma(1.25) * 1.5 * 2^0.25))^(1 / 1.5)
        noise = randn(pop_size, D)
        levy_noise = C ./ abs.(randn(pop_size, D)).^(1 / 1.5)

        PopDec = [
            elite[i] .+ noise[i, :] .* levy_noise[i, :] .* alpha .* (elite[i] .- PopDec[i] .* 2 * iter / max_iter)
            for i in 1:pop_size
        ]

        # 记录迭代点
        push!(iter_points, best)

        # 记录收敛曲线
        push!(convergence_curve, f(best))

        # 更新种群
        population = PopDec
        fitness = [f(x) for x in population]
        best_idx = argmin(fitness)
        best = copy(population[best_idx])

        verbose && println("Iter: $iter, Best fitness: $(fitness[best_idx])")
    end

    return best, iter_points, convergence_curve
end
