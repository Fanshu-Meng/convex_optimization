# Booth 函数
function booth(x)
    return (x[1] + 2 * x[2] - 7)^2 + (2 * x[1] + x[2] - 5)^2
end

# Booth 函数的梯度
function ∇booth(x)
    df_dx1 = 2 * (x[1] + 2 * x[2] - 7) + 4 * (2 * x[1] + x[2] - 5)
    df_dx2 = 4 * (x[1] + 2 * x[2] - 7) + 2 * (2 * x[1] + x[2] - 5)
    return [df_dx1, df_dx2]
end

# Booth 函数的海森矩阵
function ∇²booth(x)
    d2f_dx1dx1 = 2 + 8
    d2f_dx1dx2 = 4 + 2
    d2f_dx2dx1 = 4 + 2
    d2f_dx2dx2 = 8 + 2
    return [d2f_dx1dx1 d2f_dx1dx2; d2f_dx2dx1 d2f_dx2dx2]
end    