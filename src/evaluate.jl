""" Return number of switches. """
function n_switches(results)
    return length(results["switches"]["technology"])
end

""" Return number of reswitches. """
function n_reswitches(results)
    at = Float64[]
    # Ensure that we use technology of final switch as well
    technology = [x[1].second for x in results["switches"]["technology"]]
    push!(technology, results["switches"]["technology"][end][2].second)
    for (i, tech1) in enumerate(technology[begin:(end - 1)])
        for tech2 in technology[(i + 1):end]
            if tech1 == tech2
                push!(at, results["switches"]["technology"][i][1].first)
            end
        end
    end
    # It is possible that we reswitch multiple times with the same technology.
    return unique(at)
end

""" Return the number of cases in each of the four switching possibilities. """
function switch_cases(results)
    # TODO: check for switches in adjacent profit rates -> whether grid is fine enough
    κ_down_labour_up = κ_down_labour_down = κ_up_labour_up = κ_up_labour_down = 0
    sum_of_multiple_switches = 0
    for i in eachindex(results["switches"]["capital_intensities"])
        tech = results["switches"]["technology"][i]
        sanity_check = sum(tech[1].second .!= tech[2].second)
        (sanity_check != 1) ? sum_of_multiple_switches += 1 : nothing
        tech_index_switches = findall(tech[1].second .!= tech[2].second)
        ci = results["switches"]["capital_intensities"][i]
        chii = results["switches"]["l"][i]
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
    return sum_of_multiple_switches, κ_down_labour_up, κ_down_labour_down, κ_up_labour_up, κ_up_labour_down
end

function switch_info(results)
    # Replication number of switches
    num_switches = n_switches(results)
    @info "Number of switches: ", num_switches

    # Replication number of reswitches
    @info "Number of reswitches: ", n_reswitches(results)

    # Replication intensity cases
    found_cases = switch_cases(results)
    @info "There have been $(found_cases[1]) instances of switches that are not piecemeal."
    cases_perc = round.(found_cases ./ num_switches .* 100, digits = 2)
    @info "Capital intensity-reducing, labour-increasing: $(found_cases[2]) cases ($(cases_perc[2])%)"
    @info "Capital intensity-reducing, labour-reducing: $(found_cases[3]) cases ($(cases_perc[3])%)"
    @info "Capital intensity-increasing, labour-increasing: $(found_cases[4]) cases ($(cases_perc[4])%)"
    @info "Capital intensity-increasing, labour-reducing: $(found_cases[5]) cases ($(cases_perc[5])%)"
end

@userplot WageCurves

@recipe function f(env::WageCurves; stepsize=0.001, switches=true)
    results = env.args[1]
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
            curves_countries[i, j] = compute_w(A_curve, B_curve, d, l_curve, r) # FIXME: Accurate with big(r) (later as well)
        end
    end
    curves_env = zeros(length(grid_R), length(techs) + 1)
    r_env = zeros(size(curves_env, 2))
    A_curve = A[1:n_techs, techs[1][1].second]
    B_curve = B[1:n_techs, techs[1][1].second]
    l_curve = l[techs[1][1].second]
    r_env[1] = compute_R(real_eigvals(A_curve))
    for (i, r) in enumerate(grid_R)
        w = compute_w(A_curve, B_curve, d, l_curve, r)
        curves_env[i, 1] = r <= r_env[1] ? w : 0.0
    end
    for (j, tech) in enumerate(techs)
        A_curve = A[1:n_techs, tech[2].second]
        B_curve = B[1:n_techs, tech[2].second]
        l_curve = l[tech[2].second]
        r_env[j + 1] = compute_R(real_eigvals(A_curve))
        for (i, r) in enumerate(grid_R)
            w = compute_w(A_curve, B_curve, d, l_curve, r)
            curves_env[i, j + 1] = r <= r_env[j + 1] ? w : 0.0
        end
    end
    labels="technique " .* string.(axes(curves_env, 2))
    n_sw = n_switches(results)
    for i in axes(curves_env, 2)
        grid = 1:findlast(x -> x <= r_env[i], grid_R)
        @series begin
            title := "Wage curves with $n_sw switches"
            legend := :outertopright
            xlabel --> "r"
            ylabel --> "w"
            subplot := 1
            name --> labels[i]
            linewidth --> 3
            grid_R[grid], curves_env[grid, i]
        end
    end
    if switches
        @series begin
            title := "Wage curves with $n_sw switches"
            legend := :outertopright
            xlabel --> "r"
            ylabel --> "w"
            subplot := 1
            linewidth --> 0.5
            color --> :grey
            label --> nothing
            seriestype  := :vline
            [tech[1].first for tech in techs]
        end
    end
end

@userplot Env

@recipe function f(env::Env; stepsize=0.001, switches=true)
    results = env.args[1]
    techs = results["switches"]["technology"]
    n_techs = length(techs[end][2].second)
    A = results["A"]
    B = results["B"]
    l = vec(results["l"])
    d = results["d"][1:n_techs]
    max_R = compute_R(maximum(real_eigvals(A[1:n_techs, techs[end][2].second]')))
    grid_R = 0:stepsize:max_R
    curve = zeros(length(grid_R))
    lower_r = 0.0
    for (j, tech) in enumerate(techs)
        upper_r = techs[j][1].first
        A_curve = A[1:n_techs, tech[1].second]
        B_curve = B[1:n_techs, tech[1].second]
        l_curve = l[tech[1].second]
        for (i, r) in enumerate(grid_R)
            if lower_r <= r <= upper_r
                w = compute_w(A_curve, B_curve, d, l_curve, r)
                curve[i] = w
            end
        end
        lower_r = upper_r
    end
    for (i, r) in enumerate(grid_R)
        A_curve = A[1:n_techs, techs[end][end].second]
        B_curve = B[1:n_techs, techs[end][end].second]
        l_curve = l[techs[end][end].second]
        if techs[length(techs)][2].first <= r
            w = compute_w(A_curve, B_curve, d, l_curve, r)
            curve[i] = w
        end
    end
    n_sw = n_switches(results)
    @series begin
        title := "Wage curves with $n_sw switches"
        legend := :outertopright
        xlabel --> "r"
        ylabel --> "w"
        subplot := 1
        label --> "envelope"
        linewidth --> 3
        seriestype  := :line
        grid_R, curve
    end
end