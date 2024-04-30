""" Return total number of switches. """
function n_switches(results_pairs)
    switches = 0
    total_comparisons = 0
    for (_, val1) in results_pairs
        for (_, val2) in val1
            total_comparisons += 1
            switches += length(val2["switches"]["capital_intensities"])
        end
    end
    println("There are $total_comparisons pairwise comparisons.")
    return switches
end

""" Return number of reswitches. """
function n_reswitches(results_pairs)
    reswitches = 0
    for (country1, val1) in results_pairs
        for (country2, val2) in val1
            # Ensure that we use technology of final switch as well
            technology = [x[1].second for x in val2["switches"]["technology"]]
            push!(technology, val2["switches"]["technology"][end][2].second)
            # println(technology)
            for (i, tech1) in enumerate(technology[begin:(end - 1)])
                for tech2 in technology[(i + 1):end]
                    if tech1 == tech2
                        reswitches += 1
                        println("Reswitching in: ", "$country1", "/", "$country2",
                            ", r = $(val2["switches"]["technology"][i][1].first)")
                        break
                    end
                end
            end
        end
    end
    return reswitches
end

""" Return the number of cases in each of the four switching possibilities. """
function switch_cases(results_pairs)
    κ_down_labour_up = κ_down_labour_down = κ_up_labour_up = κ_up_labour_down = 0
    sum_of_multiple_switches = 0
    for (_, val1) in results_pairs
        for (_, val2) in val1
            for i in eachindex(val2["switches"]["capital_intensities"])
                tech = val2["switches"]["technology"][i]
                sanity_check = sum(tech[1].second .!= tech[2].second)
                (sanity_check != 1) ? sum_of_multiple_switches += 1 : nothing
                tech_index_switches = findall(tech[1].second .!= tech[2].second)
                ci = val2["switches"]["capital_intensities"][i]
                chii = val2["switches"]["l"][i]
                labour_up = any(map(i -> chii[1].second[i] <
                            chii[2].second[i], tech_index_switches))
                if ci[1].second >= ci[2].second
                    if labour_up
                        κ_down_labour_up += 1
                    else
                        κ_down_labour_down += 1
                    end
                else
                    if labour_up
                        κ_up_labour_up += 1
                    else
                        κ_up_labour_down += 1
                    end
                end
            end
        end
    end
    println("There have been $sum_of_multiple_switches instances of switches that are not piecemeal.")
    return κ_down_labour_up, κ_down_labour_down, κ_up_labour_up, κ_up_labour_down
end

function plot_wage_curves(results, stepsize=0.001)
    techs = results["switches"]["technology"]
    n_techs = length(techs[end][2].second)
    A = results["A"]
    B = results["B"]
    l = vec(results["l"])
    n_countries = div(size(A, 2), size(A, 1))
    d = results["d"][1:n_techs]
    max_R = compute_R(maximum(real_eigvals(A[1:n_techs, techs[end][2].second]')))
    grid_R = 0:stepsize:max_R
    curves_countries = zeros(length(grid_R), n_countries)
    for j in 1:n_countries
        tech = (1:n_techs) .+ (j - 1) * size(A, 1)
        A_curve = A[1:n_techs, tech]
        B_curve = B[1:n_techs, tech]
        l_curve = l[tech]
        for (i, r) in enumerate(grid_R)
            curves_countries[i, j] = compute_w(A_curve, B_curve, d, l_curve, r)
        end
    end
    curves_env = zeros(length(grid_R), length(techs) + 1)
    r_env = zeros(size(curves_env, 2))
    A_curve = A[1:n_techs, techs[1][1].second]
    B_curve = B[1:n_techs, techs[1][1].second]
    l_curve = l[techs[1][1].second]
    r_env[1] = compute_R(maximum(real_eigvals(A_curve)))
    for (i, r) in enumerate(grid_R)
        w = compute_w(A_curve, B_curve, d, l_curve, r)
        curves_env[i, 1] = r <= r_env[1] ? w : 0.0
    end
    for (j, tech) in enumerate(techs)
        A_curve = A[1:n_techs, tech[2].second]
        B_curve = B[1:n_techs, tech[2].second]
        l_curve = l[tech[2].second]
        r_env[j + 1] = compute_R(maximum(real_eigvals(A_curve)))
        for (i, r) in enumerate(grid_R)
            w = compute_w(A_curve, B_curve, d, l_curve, r)
            curves_env[i, j + 1] = r <= r_env[j + 1] ? w : 0.0
        end
    end
    p = plot(legend=:outertopright)
    labels="tech " .* string.(axes(curves_env, 2))
    for i in axes(curves_env, 2)
        grid = 1:findlast(x -> x <= r_env[i], grid_R)
        plot!(grid_R[grid], curves_env[grid, i], name=labels[i], xlabel="r", ylabel="w", lw = 3)
    end
    # p = plot!(grid_R, curves_env, label="country " .* string.(axes(curves_countries, 2)), color = :gray, lw = 2)
    vline!([tech[1].first for tech in techs], label=nothing, color=:grey)
    return p
end