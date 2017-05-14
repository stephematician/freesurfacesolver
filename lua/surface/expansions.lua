local d_arclength_LIP_estimate
local first_order_stepsize
local second_order_stepsize

d_arclength_LIP_estimate = function(t, y)

    local ds = 0
    local d2s = 0

    for k = 1, #y do
        local f, df, d2f

        f = {y[k][1],
             y[k][2],
             y[k][3]}
        df = -f[1] * t[2] / (t[1] * (t[1] - t[2])) +
             -f[2] * t[1] / (t[2] * (t[2] - t[1])) +
             -f[3] * (t[1] + t[2]) / (t[1] * t[2])
        d2f = 2 * f[1] / (t[1] * (t[1] - t[2]) ) +
              2 * f[2] / (t[2] * (t[2] - t[1])) +
              2 * f[3] / (t[1] * t[2])

        ds = ds + df * df
        d2s = d2s + (df * d2f * d2f)

    end

    ds = math.sqrt(ds)
    d2s = d2s / ds

    return ds, d2s

end

first_order_stepsize = function(step_history, curve_history, ds)

    local ds1 = 0

    for k = 1, #curve_history.peek(1) do
        ds1 = ds1 + math.pow(curve_history.peek(1)[k] -
                                 curve_history.peek(2)[k], 2)
    end
    ds1 = math.sqrt(ds1)

    return step_history.peek(1) * ds / ds1

end

second_order_stepsize = function(step_history, curve_history, ds)
    
    local y = {}
    local t = {-(step_history.peek(2) + step_history.peek(1)), 
               -(step_history.peek(1)),
               0}

    for k = 1, #surf.length_vars do
        y[k][1] = curve_history.peek(3)[k]
        y[k][2] = curve_history.peek(2)[k]
        y[k][3] = curve_history.peek(1)[k]
    end

    local LIP_ds, LIP_d2s = d_arclength_LIP_estimate(t, y)
    local a = LIP_d2s / 2
    local b = LIP_ds
    local c = -ds

    local step = -1/0

    if b*b - 4*a*c > 0 then
        step = ((-b + math.sqrt(b*b - 4*a*c))/(2*a))
    end

    return step

end

return {
    first_order_stepsize = first_order_stepsize,
    second_order_stepsize = second_order_stepsize
}

