using Polynomials
using Plots
using LaTeXStrings

function getv(n, l)
    if n > l
        a = zeros(n+1)
        a[1] = 1
        for j = 1:n
            a[j+1] = (2*(j+l) - 2*n)/(j*(j+2*l+1))*a[j]
        end
        return Polynomial(a)
    else
        return Polynomial([0])
    end
end

function getu(n, l)
    v = getv(n, l)
    u(ρ) = ρ^(l+1)*exp(-ρ)*v(ρ)
    return u
end

function getR(n, l)
    u = getu(n, l)
    R(r) = u(r/n)/r
    return R
end

function getlegendre(l)
    a = zeros(l+1)

    start = 1
    if l % 2 != 0
        start = 2
    end
    a[start] = 1

    β = l*(l+1)
    for j in (start+2):2:(l+1)
        a[j] = ((j-3)*(j-2) - β)/((j-1)*(j-2))*a[j-2]
    end
    
    unnormalized = Polynomial(a)
    return unnormalized/unnormalized(1)
end

function getassoclegendre(l, m)
    P_l = getlegendre(l)
    P_lm(x) = (1-x^2)^(abs(m)/2) * derivative(P_l, abs(m))(x)
    return P_lm
end

function getΘ(l, m)
    P_lm = getassoclegendre(l, m)
    Θ(θ) = P_lm(cos(θ))
    return Θ
end

function getΦ(m)
    if m >= 0
        Φ(ϕ) = cos(m*ϕ)
        return Φ
    else
        Φ(ϕ) = sin(m*ϕ)
        return Φ
    end
end

function getanalyticalψ(n, l, m; L=100.0)
    R = getR(n, l)
    Θ = getΘ(l, m)
    Φ = getΦ(m)

    rnorm = solve(IntegralProblem((r, p) -> abs(R(r))^2*r^2, 0.0, Inf), HCubatureJL(); reltol = 1e-8, abstol = 1e-8).u
    thetanorm = solve(IntegralProblem((r, p) -> abs(Θ(r))^2*sin(r), 0.0, float(pi)), HCubatureJL(); reltol = 1e-8, abstol = 1e-8).u
    phinorm = solve(IntegralProblem((r, p) -> abs(Φ(r))^2, 0.0, 2*float(pi)), HCubatureJL(); reltol = 1e-8, abstol = 1e-8).u

    wf = [R(r)*Θ(θ)*Φ(ϕ) for r in range(0, L, 11)[2:end], θ in range(0, pi, 11)[1:end-1] .+ 0.5*(pi/10), ϕ in range(0, 2*pi, 11)[1:end-1]]
    flip = sign(wf[findmax(abs.(wf))[2]])
    wfnorm_analytical = sqrt(rnorm*thetanorm*phinorm)*flip

    function ψ(r, θ, ϕ)
        return R(r) * Θ(θ) * Φ(ϕ) / wfnorm_analytical
    end
    return ψ
end