using Plots
gr()  # 使用 GR 后端

# 包含算法和测试函数文件
include("algorithms/simplex.jl")
include("algorithms/gradient_descent.jl")
include("algorithms/subgradient_descent.jl")
include("algorithms/conjugate_direction.jl")
include("algorithms/conjugate_gradient.jl")
include("algorithms/quasi_newton.jl")
include("algorithms/momentum.jl")
include("algorithms/newton.jl")
include("algorithms/damped_newton.jl")
include("algorithms/cubic_regularized_newton.jl")
include("algorithms/admm.jl")
include("algorithms/krylov.jl")

include("algorithms/DOA.jl")

include("test_functions/rosenbrock.jl")
include("test_functions/booth.jl")

# 定义初始点
x0 = [-1.0, -1.0]

# 定义优化算法列表
algorithms = [
    ("simplex", simplex),
    ("gradient_descent", gradient_descent),
    ("subgradient_descent", subgradient_descent),
    ("conjugate_direction", conjugate_direction),
    ("conjugate_gradient", conjugate_gradient),
    ("quasi_newton", quasi_newton),
    ("momentum", momentum),
    ("newton", newton),
    ("damped_newton", damped_newton),
    ("cubic_regularized_newton", cubic_regularized_newton),
    ("admm", admm),
    ("krylov", krylov),
    ("DOA", DOA)
]

# 定义测试函数信息
test_functions = [
    (rosenbrock, ∇rosenbrock, ∇²rosenbrock, "rosenbrock", [1.0, 1.0]),
    (booth, ∇booth, ∇²booth, "booth", [1.0, 3.0])
]

# 创建保存目录
if !isdir("figures")
    mkdir("figures")
end

txt_path = joinpath("figures", "log.txt")

# 遍历算法
for (alg_name, alg_func) in algorithms
    println("Running algorithm: $alg_name")
    alg_folder = joinpath("figures", alg_name)
    if !isdir(alg_folder)
        mkdir(alg_folder)
    end

    for (f, ∇f, ∇²f, func_name, correct_answer) in test_functions
        func_folder = joinpath(alg_folder, func_name)
        if !isdir(func_folder)
            mkdir(func_folder)
        end

        func = (f, ∇f, ∇²f) -> alg_func(f, ∇f, ∇²f, x0)
        start_time = time()
        result, iter_points, convergence_curve = func(f, ∇f, ∇²f)
        run_time = time() - start_time

        function_value = f(result)
        error = norm(result - correct_answer)
        iterations = length(iter_points)

        # 提取坐标
        xs = [p[1] for p in iter_points]
        ys = [p[2] for p in iter_points]
        x_opt = correct_answer

        # 绘制等高线 + 轨迹图
        plt_iter = contour(
            -1.5:0.1:2, -1.5:0.1:3.5, 
            (x1, x2) -> f([x1, x2]), 
            levels = 50, 
            fill = true, 
            colormap = :cividis, 
            linecolor = :black, 
            linewidth = 0.5, 
            title = "Iteration Path",
            xlabel = "x1", 
            ylabel = "x2",
            legend = :bottomright
        )

        # 虚线轨迹
        plot!(plt_iter, xs, ys, linestyle = :dash, color = :red, label = "$alg_name path")

        # 起点
        scatter!(plt_iter, [xs[1]], [ys[1]], color = :blue, marker = (:star5, 10), label = "Start")

        # 终点
        scatter!(plt_iter, [xs[end]], [ys[end]], color = :green, marker = (:circle, 6), label = "End")

        # 最优点
        scatter!(plt_iter, [x_opt[1]], [x_opt[2]], color = :black, marker = (:diamond, 8), label = "Optimal Point")

        # 保存轨迹图
        savefig(plt_iter, joinpath(func_folder, "iterative_curve.png"))

        # 收敛曲线图
        plt_conv = plot(
            convergence_curve,
            title = "Convergence Curve",
            xlabel = "Iteration", 
            ylabel = "Function Value",
            legend = false
        )
        savefig(plt_conv, joinpath(func_folder, "convergence_curve.png"))

        # 写入日志
        open(txt_path, "a+") do f
            println(f, "Algorithm: $alg_name")
            println(f, "Test Function: $func_name")
            println(f, "Optimal Solution: $(join(result, ", "))")
            println(f, "Function Value: $(round(function_value, digits=4))")
            println(f, "Error (L2 Norm): $(round(error, digits=4))")
            println(f, "Iterations: $iterations")
            println(f, "Run Time: $(round(run_time, digits=4)) seconds")
            println(f, "\n")
        end
    end
end
