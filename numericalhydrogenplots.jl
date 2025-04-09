using Plots
using LaTeXStrings
using LinearAlgebra
using SparseArrays
using Polynomials
using Integrals

include("analyticalhydrogenplots.jl")

function generateΦmatrix(N, h)
    offdiag = fill(1/h^2, N-1)
    diagonal = fill(-2/h^2, N)
    A = Matrix(SymTridiagonal(diagonal, offdiag))
    A[1, N] = 1/h^2
    A[N, 1] = 1/h^2
    return A
end

function solveΦ(m, N)
    h = 2*pi/(N+1)

    A = generateΦmatrix(N+1, h)
    eigendecomp = eigen(A)
    sol1 = sign(eigendecomp.vectors[:, end-2*abs(m)][1]) * eigendecomp.vectors[:, end-2*abs(m)]             # need to pick sin or cos. m < 0 => sin, m > 0 => cos
    sol2 = sign(eigendecomp.vectors[:, end-2*abs(m)+(m != 0 ? 1 : 0)][1]) * eigendecomp.vectors[:, end-2*abs(m)+(m != 0 ? 1 : 0)]
    solutions = [max(sol1, sol2), min(sol1, sol2)]
    Φpoints = solutions[m > 0 ? 1 : 2]
    function Φ(ϕ)
        i = Int(ϕ ÷ h)
        s = ϕ/h - i
        if i == 0
            return Φpoints[end]
        end
        return ((1-s) * Φpoints[i] + s * Φpoints[i+1])
    end

    eigenvalue = eigendecomp.values[end-2*abs(m)+(m>0 ? 1 : 0)]
    return (Φ, eigenvalue)
end

function generateΘmatrix(θs, m, h)
    subdiag = @. 1/h^2 - cot(θs[2:end])/(2*h)
    diagonal = @. -2/h^2 - m^2/sin(θs)^2
    superdiag = @. 1/h^2 + cot(θs[1:end-1])/(2*h)

    A = Matrix(Tridiagonal(subdiag, diagonal, superdiag))

    A[1, 1] = A[1, 1] + 1/h^2 - cot(θs[1])/(2*h)
    A[length(θs), length(θs)] = A[length(θs), length(θs)] + 1/h^2 + cot(θs[end])/(2*h)

    #A[1, length(θs)] = 1/h^2 - cot(θs[1])/(2*h)
    #A[length(θs), 1] = 1/h^2 + cot(θs[end])/(2*h)
    return A
end

function solveΘ(l, m, N)
    h = pi/(N+1)

    θs = range(0, pi, N+2)[1:end-1] .+ 0.5*h
    A = generateΘmatrix(θs, abs(m), h)

    eigendecomp = eigen(A)
    Θpoints = eigendecomp.vectors[:, end-l+abs(m)]
    eigenvalue = eigendecomp.values[end-l+abs(m)]

    function Θ(θ)
        i = Int((θ-h/2) ÷ h)
        s = (θ-h/2)/h - i

        return ((1-s) * Θpoints[i > N ? 2*N-i+2 : i+1] + s * Θpoints[i+1 > N ? 2*N-i+1 : i+2])
    end
    return (Θ, eigenvalue)
end

function generateRmatrix(rs, l, h)
    subdiag = @. -1/(2*h^2) + 1/(2*rs[2:end]*h)
    diagonal = @. 1/h^2 + l*(l+1)/(2*rs^2) - 1/rs
    superdiag = @. -1/(2*h^2) - 1/(2*rs[1:end-1]*h)

    A = Tridiagonal(subdiag, diagonal, superdiag)

    return A
end

function solveR(n, l, N, L)
    h = L/(N+1)

    rs = range(0, L, N+2)[2:end-1]
    A = generateRmatrix(rs, l, h)

    eigendecomp = eigen(A)
    Rpoints = [0; eigendecomp.vectors[:, n-l]; 0]
    eigenvalue = eigendecomp.values[n-l]

    function R(r)
        i = Int(r ÷ h)
        s = r/h - i
        if i == N+1
            return Rpoints[end]
        end
        return ((1-s) * Rpoints[i+1] + s * Rpoints[i+2])
    end
    return (R, eigenvalue)
end

function solveψ(n, l, m, N, L)
    (Φ, Φeigenvalue) = solveΦ(m, N)
    (Θ, Θeigenvalue) = solveΘ(l, m, N)
    (R, Reigenvalue) = solveR(n, l, N, L)

    rnorm = solve(IntegralProblem((r, p) -> abs(R(r))^2*r^2, 0.0, L), HCubatureJL(); reltol = 1e-8, abstol = 1e-8).u
    thetanorm = solve(IntegralProblem((r, p) -> abs(Θ(r))^2*sin(r), 0.0, float(pi)), HCubatureJL(); reltol = 1e-8, abstol = 1e-8).u
    phinorm = solve(IntegralProblem((r, p) -> abs(Φ(r))^2, 0.0, 2*float(pi)), HCubatureJL(); reltol = 1e-8, abstol = 1e-8).u

    wf = [R(r)*Θ(θ)*Φ(ϕ) for r in range(0, L, 11)[2:end], θ in range(0, pi, 11)[1:end-1] .+ 0.5*(pi/10), ϕ in range(0, 2*pi, 11)[1:end-1]]
    flip = sign(wf[findmax(abs.(wf))[2]])
    wfnorm = sqrt(rnorm*thetanorm*phinorm)*flip

    function ψ(r, θ, ϕ)
        return R(r) * Θ(θ) * Φ(ϕ) / wfnorm
    end

    numericalm = sqrt(abs(Φeigenvalue))
    numericall = (-1 + sqrt(4*abs(Θeigenvalue) + 1))/2
    numericalE = Reigenvalue
    numericaln = sqrt(abs(1/(2*numericalE)))
    return (ψ, numericalE, numericaln, numericall, numericalm)
end

function solve_hydrogen(n, l, m, N, L; plotting=true, plotL=10)
    h_r = L/(N+1)
    h_theta = pi/(N+1)

    (ψ, numericalE, numericaln, numericall, numericalm) = solveψ(n, l, m, N, L)
    println("numerical:  E = ", numericalE, "; n = ", numericaln, "; l = ", numericall, "; m = ", numericalm)
    println("analytical: E = ", -0.5/n^2, "; n = ", n, "; l = ", l, "; m = ", m)

    analyticalψ = getanalyticalψ(n, l, m)

    stride = Int(N/min(N, 512))
    error_phis = (range(0, 2*pi, N+1)[1:end-1])[1:stride:end]
    error_thetas = (range(0, pi, N+2)[1:end-1] .+ 0.5*h_theta)[1:stride:end]
    error_rs = (range(0, L, N+2)[2:end-1])[1:stride:end]

    max_abs_error = 0
    max_rel_error = 0
    for r in error_rs, theta in error_thetas, phi in error_phis
        analytical = analyticalψ(r, theta, phi)
        numerical = ψ(r, theta, phi)

        abs_error = (analytical-numerical)
        rel_error = abs_error/(analytical >= 1e-12 ? analytical : Inf)

        max_abs_error = max(max_abs_error, abs(abs_error))
        max_rel_error = max(max_rel_error, abs(rel_error))
    end

    println("maximum absolute error: ", max_abs_error)
    println("maximum relative error: ", max_rel_error)

    l2error = solve(IntegralProblem((r, p) -> abs(ψ(r[1], r[2], r[3]) - analyticalψ(r[1], r[2], r[3]))^2*r[1]^2*sin(r[2]), [h_r, 0.0, 0.0], [L, float(pi), 2*float(pi)]), HCubatureJL(); reltol = 1e-6, abstol = 1e-6).u
    println("L^2 error: ", l2error)

    if plotting
        plotN = min(N, 128)
        stride = Int(N/plotN)
        plot_phis = [fill(0, plotN+1); fill(pi, plotN+1)]
        plot_thetas = (range(0, pi, N+2)[1:end-1] .+ 0.5*h_theta)[1:stride:end]
        plot_thetas = [plot_thetas; plot_thetas[end:-1:1]]
        rstride = Int(floor(N*(plotL/L)/min(N*(plotL/L), 128)))
        plot_rs = (range(0, L, N+2)[2:end-1])[1:rstride:end]
        plot_rs = plot_rs[plot_rs .<= plotL]

        analyticalwavefunction = [analyticalψ(r, angles[1], angles[2]) for r in plot_rs, angles in zip(plot_thetas, plot_phis)]

        wavefunction = [ψ(r, angles[1], angles[2]) for r in plot_rs, angles in zip(plot_thetas, plot_phis)]
        probs = @. abs(wavefunction)^2

        abs_error = wavefunction-analyticalwavefunction
        abs_errorrange = maximum(abs, abs_error)

        orb = heatmap(range(0, 2*pi, 2*plotN+2), plot_rs, probs*1e3, projection=:polar, grid=false, axis=false, color=cgrad(:inferno, scale=:exp), rightmargin=10Plots.mm, dpi=1000,xtickfontsize=14,ytickfontsize=14,xguidefontsize=14,yguidefontsize=14,legendfontsize=14, leftmargin=-30Plots.mm)
        annotate!(1.7, -1.05, L"\times 10^{-3}")
        scatter!([0], [0], label="nucleus", mc=:green, m=:cross, ms=6, markerstrokewidth=2, legend=:topright)
        savefig("plots/numerical/numericalprob_$(n)$(l)$(m).png")

        waverange = maximum(abs, wavefunction)
        wave = heatmap(range(0, 2*pi, 2*plotN+2), plot_rs, wavefunction*1e2, projection=:polar, grid=false, axis=false, color=:seaborn_icefire_gradient, clim=(-waverange*1e2, waverange*1e2), rightmargin=10Plots.mm, dpi=1000,xtickfontsize=14,ytickfontsize=14,xguidefontsize=14,yguidefontsize=14,legendfontsize=14, leftmargin=-30Plots.mm)
        annotate!(1.7, -1.05, L"\times 10^{-2}")
        scatter!([0], [0], label="nucleus", mc=:green, m=:cross, ms=6, markerstrokewidth=2, legend=:topright)
        savefig("plots/numerical/numericalwave_$(n)$(l)$(m).png")

        err = heatmap(range(0, 2*pi, 2*plotN+2), plot_rs, abs_error*1e5, projection=:polar, grid=false, axis=false, color=:seaborn_icefire_gradient, clim=(-abs_errorrange*1e5, abs_errorrange*1e5), rightmargin=10Plots.mm, dpi=1000,xtickfontsize=14,ytickfontsize=14,xguidefontsize=14,yguidefontsize=14,legendfontsize=14, leftmargin=-30Plots.mm)
        annotate!(1.7, -1.05, L"\times 10^{-5}")
        scatter!([0], [0], label="nucleus", mc=:green, m=:cross, ms=6, markerstrokewidth=2, legend=:topright)
        savefig("plots/numerical/numericalabserror_$(n)$(l)$(m).png")
    end
end

gr()
n, l, m = 2, 1, 0
L = 100.0

N = 2048
solve_hydrogen(n, l, m, N, L, plotting=true, plotL=15)

for i = 6:12
    N = 2^i
    println()
    println("N = ", N, ":")
    solve_hydrogen(n, l, m, N, L, plotting=false)
end