module Projection
    using JuMP
    using Random
    using Manifolds

    function get_model(solver)
        model = direct_model(solver.Optimizer())
        set_silent(model)
        return model
    end

    function bounds(T::Union{Matrix{Float64}, Matrix{Int64}}, t::Union{Vector{Float64}, Vector{Int64}}, Q::Model)
        upper_bounds = []
        lower_bounds = []
        p = size(T, 1)
        q = size(T, 2)
        y = all_variables(Q)
        for i in 1:p
            @objective(Q, Max, T[i,1:q]' * y)
            optimize!(Q)
            push!(upper_bounds, objective_value(Q) + t[i])

            @objective(Q, Min, T[i,1:q]' * y)
            optimize!(Q)
            push!(lower_bounds, objective_value(Q) + t[i])
        end
        return lower_bounds, upper_bounds
    end

    function get_separator(A::Union{Matrix{Float64}, Matrix{Int64}}, b::Union{Vector{Float64}, Vector{Int64}}, T::Union{Matrix{Float64}, Matrix{Int64}}, t::Union{Vector{Float64}, Vector{Int64}}, solver::Module)
        p = size(T,1)
        sep = get_model(solver)
        @variable(sep,  a[1:p])
        @variable(sep, -1 <= β <= 1)
        @variable(sep, 0 <= λ[1:size(A,1)])
        @constraint(sep, λ' * A .== a' * T)
        @constraint(sep, λ' * b + a' * t <= β)
        return sep
    end

    function get_cut(sep::Model, tilde_x::Vector{Float64}, p::Int64)
        sep_vars = all_variables(sep)
        a = sep_vars[1:p]
        β = sep_vars[p+1]
        @objective(sep, Max, tilde_x' * a - β)
        optimize!(sep)
        sep_solution = value.(all_variables(sep))
        return sep_solution[1:p], sep_solution[p+1]
    end

    function approximate_projection(A::Union{Matrix{Float64}, Matrix{Int64}}, b::Union{Vector{Float64}, Vector{Int64}}, T::Union{Matrix{Float64}, Matrix{Int64}}, t::Union{Vector{Float64}, Vector{Int64}}, threshold::Float64, termcount::Integer, solver::Module)
        totalIter = 0
        gapLarge = true
        constrList = []         
        p = size(T, 1)
        q = size(T, 2)

        Q = get_model(solver)
        @variable(Q, y[1:q])
        @constraint(Q, A * y .<= b)

        P = get_model(solver)
        @variable(P, x[1:p])

        sep = get_separator(A, b, T, t, solver)

        lower_bounds, upper_bounds = bounds(T, t, Q)
        for i in 1:p
            lhs = zeros(p)
            lhs[i] = 1
            rhs = upper_bounds[i]
            push!(constrList, (lhs, rhs))
            @constraint(P, lhs' * x <= rhs)

            lhs = zeros(p)
            lhs[i] = -1
            rhs = -lower_bounds[i]
            push!(constrList, (lhs, rhs))
            @constraint(P, lhs' * x <= rhs)
        end

        sphere = Sphere(p - 1)
        sample_space = uniform_distribution(sphere, zeros(representation_size(sphere)))
        while gapLarge
            gapLarge = false
            for i in 1:termcount
                totalIter = totalIter + 1
                rand_obj = broadcast(abs, rand(sample_space, 1)[1])
                @objective(P, Max, rand_obj' * x)
                optimize!(P)
                @objective(Q, Max, rand_obj' * (T * y + t))
                optimize!(Q)
                if (objective_value(P) - objective_value(Q)) / objective_value(Q) > threshold
                    gapLarge = true
                    tilde_x = value.(all_variables(P))
                    lhs, rhs = get_cut(sep, tilde_x, p)
                    @constraint(P, lhs' * x <= rhs)
                    push!(constrList, (lhs, rhs))
                end
            end
        end
        return constrList
    end
end