# Rosenbrock Banana 函数
function rosenbrock(x; a = 1, b = 5)
    return (a - x[1])^2 + b * (x[2] - x[1]^2)^2
end

# Rosenbrock Banana 函数的梯度
function ∇rosenbrock(x; a = 1, b = 5)
    df_dx1 = -2 * (a - x[1]) - 4 * b * x[1] * (x[2] - x[1]^2)
    df_dx2 = 2 * b * (x[2] - x[1]^2)
    return [df_dx1, df_dx2]
end

# Rosenbrock Banana 函数的海森矩阵
function ∇²rosenbrock(x; a = 1, b = 5)
    d2f_dx1dx1 = 2 - 4 * b * (x[2] - 3 * x[1]^2)
    d2f_dx1dx2 = -4 * b * x[1]
    d2f_dx2dx1 = -4 * b * x[1]
    d2f_dx2dx2 = 2 * b
    return [d2f_dx1dx1 d2f_dx1dx2; d2f_dx2dx1 d2f_dx2dx2]
end    