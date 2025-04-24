# 多维梯度下降法
function gradient_descent(f, ∇f, ∇²f, x0; α=0.01, max_iter=1000, tol=1e-6)
    x = x0
    iter_points = [x0]
    convergence_curve = [f(x0)]

    for iter in 1:max_iter
        g = ∇f(x)  # 计算梯度
        x_new = x - α * g  # 梯度下降更新

        push!(iter_points, x_new)
        push!(convergence_curve, f(x_new))

        if norm(x_new - x) < tol
            break
        end

        x = x_new
    end

    return x, iter_points, convergence_curve
end
