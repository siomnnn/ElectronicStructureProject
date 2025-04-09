using LinearAlgebra
using SpecialFunctions
using Plots
using Integrals
using Cubature

function SCF_iteration(N, P_0, H_core, generate_G, X, tol)
    P = P_0
    P_old = fill(Inf, size(P))
    
    K = size(P)[1]
    C = zeros(size(P))
    F = zeros(size(P))
    eigendecomp = []

    i=0
    while norm(P-P_old) >= tol
        i = i+1
        println(i, ": ", norm(P-P_old))
        P_old = P

        G = generate_G(P_old)

        F = H_core + G
        F_prime = X' * F * X

        eigendecomp = eigen(F_prime)
        C = X * eigendecomp.vectors

        P_(i, j) = 2*sum([C[i, a]*C[j, a] for a in 1:Int(N/2)])
        P = [P_(i, j) for i in 1:K, j in 1:K]
        println(eigendecomp.values)
        ground_state_energy = 0.5*sum([P_old[i, j]*(H_core[j, i] + F[j, i]) for i in 1:K, j in 1:K])
        println(ground_state_energy)
    end

    ground_state_energy = 0.5*sum([P_old[i, j]*(H_core[j, i] + F[j, i]) for i in 1:K, j in 1:K])

    return C, eigendecomp.values, ground_state_energy
end

function generate_G_from_ints(P, two_electron_integrals)
    K = size(P)[1]
    return [sum([P[λ, σ] * (two_electron_integrals[μ, ν, σ, λ] - 0.5*two_electron_integrals[μ, λ, σ, ν]) for λ in 1:K, σ in 1:K]) for μ in 1:K, ν in 1:K]
end

function generate_overlap_matrix(rs, αs)
    K = length(αs)
    overlap(i, j) = (pi/(αs[i]+αs[j]))^(3/2) * exp(-αs[i]*αs[j]*norm(rs[i, :]-rs[j, :])^2/(αs[i]+αs[j]))
    S = [(4*αs[i]*αs[j])^(3/4)/pi^(3/2)*overlap(i, j) for i in 1:K, j in 1:K]
    return S
end
function generate_H_core(rs, αs, Rs, Zs)
    K = length(αs)
    N_n = length(Zs)
    kinetic(i, j) = αs[i]*αs[j]/(αs[i] + αs[j])*(3 - 2*αs[i]*αs[j]*norm(rs[i, :]-rs[j, :])^2/(αs[i] + αs[j])) * (pi/(αs[i]+αs[j]))^(3/2) * exp(-αs[i]*αs[j]*norm(rs[i, :]-rs[j, :])^2/(αs[i]+αs[j]))
    
    F_0(t) = t == 0 ? 1 : 0.5*(pi/t)^(1/2) * erf(t^(1/2))
    nuclear(i, j) = sum([-2*pi*Zs[k]/(αs[i]+αs[j]) * exp(-αs[i]*αs[j]/(αs[i] + αs[j])*norm(rs[i, :]-rs[j, :])^2) * F_0((αs[i]+αs[j])*norm((αs[i]*rs[i, :]+αs[j]*rs[j, :])/(αs[i]+αs[j])-Rs[k, :])^2) for k in 1:N_n])
    H_core = [(4*αs[i]*αs[j])^(3/4)/pi^(3/2)*(kinetic(i, j) + nuclear(i, j)) for i in 1:K, j in 1:K]
    return H_core
end
function generate_two_electron_integrals(rs, αs)
    K = length(αs)
    F_0(t) = t == 0 ? 1 : 0.5*(pi/t)^(1/2) * erf(t^(1/2))
    two_electron(i, j, k, l) = 2*pi^(5/2)/((αs[i]+αs[j])*(αs[k]+αs[l])*(αs[i]+αs[j]+αs[k]+αs[l])^(1/2)) * 
        exp(-αs[i]*αs[j]/(αs[i]+αs[j])*norm(rs[i, :]-rs[j, :])^2-αs[k]*αs[l]/(αs[k]+αs[l])*norm(rs[k, :]-rs[l, :])^2) * 
        F_0((αs[i]+αs[j])*(αs[k]+αs[l])/(αs[i]+αs[j]+αs[k]+αs[l])*norm((αs[i]*rs[i, :]+αs[j]*rs[j, :])/(αs[i]+αs[j]) - 
        (αs[k]*rs[k, :]+αs[l]*rs[l, :])/(αs[k]+αs[l]))^2)
    two_electron_ints = [(16*αs[i]*αs[j]*αs[k]*αs[l])^(3/4)/pi^3*two_electron(i, j, k, l) for i in 1:K, j in 1:K, k in 1:K, l in 1:K]
    return two_electron_ints
end
function nuclear_energy(Rs, Zs)
    N_n = size(Rs)[1]
    energy(i, j) = (i == j ? 0 : Zs[i]*Zs[j]/norm(Rs[i, :] - Rs[j, :]))
    return sum([0.5*energy(i, j) for i in 1:N_n, j in 1:N_n])
end
function get_orbitals(C, αs, rs)
    K = length(αs)
    basis = [x -> (2*αs[i])^(3/4)/pi^(3/4)*exp(-αs[i]*norm(x - rs[i, :])^2) for i in 1:K]
    orbitals = [(x, y, z) -> sum([C[j, i]*basis[j]([x, y, z]) for j in 1:K]) for i in 1:K]
    return orbitals
end

# He
#N = 2
#αs = [0.31, 1.16, 6.36]
#K = length(αs)
#rs = zeros(K, 3)
#Rs = zeros(1, 3)
#Zs = [2]

# Be
#N = 4
#αs = [30.16787069, 5.495115306, 1.487192653, 1.314833110, 0.3055389383, 0.09937074560]
#K = length(αs)
#rs = zeros(K, 3)
#Rs = zeros(1, 3)
#Zs = [4]

# H_2 3 gaussians
#bond_length = 1.4
#N = 2
#αs = [0.168856, 0.623913, 3.42525, 0.168856, 0.623913, 3.42525]
#K = length(αs)
#rs = [bond_length/2 0 0; bond_length/2 0 0; bond_length/2 0 0; -bond_length/2 0 0; -bond_length/2 0 0; -bond_length/2 0 0]
#Rs = [bond_length/2 0 0; -bond_length/2 0 0]
#Zs = [1, 1]

# H_2 3 gaussians but twice
#bond_length = 1.4
#N = 4
#αs = [0.168856, 0.623913, 3.42525, 0.168856, 0.623913, 3.42525, 0.168856, 0.623913, 3.42525, 0.168856, 0.623913, 3.42525]
#K = length(αs)
#rs = [bond_length/2 3 0; bond_length/2 3 0; bond_length/2 3 0; -bond_length/2 3 0; -bond_length/2 3 0; -bond_length/2 3 0; bond_length/2 -3 0; bond_length/2 -3 0; bond_length/2 -3 0; -bond_length/2 -3 0; -bond_length/2 -3 0; -bond_length/2 -3 0]
#Rs = [bond_length/2 3 0; -bond_length/2 3 0; bond_length/2 -3 0; -bond_length/2 -3 0]
#Zs = [1, 1, 1, 1]

# H_2 6 gaussians
#bond_length = 1.4
#N = 2
#αs = [35.52322122, 6.513143725, 1.822142904, 0.6259552659, 0.2430767471, 0.1001124280, 35.52322122, 6.513143725, 1.822142904, 0.6259552659, 0.2430767471, 0.1001124280]
#K = length(αs)
#rs = [bond_length/2 0 0; bond_length/2 0 0; bond_length/2 0 0; bond_length/2 0 0; bond_length/2 0 0; bond_length/2 0 0; -bond_length/2 0 0; -bond_length/2 0 0; -bond_length/2 0 0; -bond_length/2 0 0; -bond_length/2 0 0; -bond_length/2 0 0]
#Rs = [bond_length/2 0 0; -bond_length/2 0 0]
#Zs = [1, 1]

# LiH
#bond_length = 1.65
#N = 4
#αs = [16.1195, 2.93620, 0.794650, 0.636289, 0.147860, 0.0480886, 0.168856, 0.623913, 3.42525]
#K = length(αs)
#rs = [bond_length/2 0 0; bond_length/2 0 0; bond_length/2 0 0; bond_length/2 0 0; bond_length/2 0 0; bond_length/2 0 0; -bond_length/2 0 0; -bond_length/2 0 0; -bond_length/2 0 0]
#Rs = [bond_length/2 0 0; -bond_length/2 0 0]
#Zs = [3, 1]

# Li2
bond_length = 5.22
N = 6
αs = [16.1195, 2.93620, 0.794650, 0.636289, 0.147860, 0.0480886, 16.1195, 2.93620, 0.794650, 0.636289, 0.147860, 0.0480886]
K = length(αs)
rs = [bond_length/2 0 0; bond_length/2 0 0; bond_length/2 0 0; bond_length/2 0 0; bond_length/2 0 0; bond_length/2 0 0; -bond_length/2 0 0; -bond_length/2 0 0; -bond_length/2 0 0; -bond_length/2 0 0; -bond_length/2 0 0; -bond_length/2 0 0]
Rs = [bond_length/2 0 0; -bond_length/2 0 0]
Zs = [3, 3]

# HeH+ 3 gaussians
#bond_length = 1.45
#N = 2
#αs = [0.168856, 0.623913, 3.42525, 0.31, 1.16, 6.36]
#K = length(αs)
#rs = [bond_length/2 0 0; bond_length/2 0 0; bond_length/2 0 0; -bond_length/2 0 0; -bond_length/2 0 0; -bond_length/2 0 0]
#Rs = [bond_length/2 0 0; -bond_length/2 0 0]
#Zs = [1, 2]

# H3+ 3 gaussians
#radius = 1.65/sqrt(3)
#N = 2
#αs = [0.168856, 0.623913, 3.42525, 0.168856, 0.623913, 3.42525, 0.168856, 0.623913, 3.42525]
#K = length(αs)
#rs = [radius 0 0; radius 0 0; radius 0 0; -radius/2 sqrt(3)/2*radius 0; -radius/2 sqrt(3)/2*radius 0; -radius/2 sqrt(3)/2*radius 0; -radius/2 -sqrt(3)/2*radius 0; -radius/2 -sqrt(3)/2*radius 0; -radius/2 -sqrt(3)/2*radius 0]
#Rs = [radius 0 0; -radius/2 sqrt(3)/2*radius 0; -radius/2 -sqrt(3)/2*radius 0]
#Zs = [1, 1, 1]

# BeH2 3 gaussians
#N = 6
#αs = [30.16787069, 5.495115306, 1.487192653, 1.314833110, 0.3055389383, 0.09937074560, 0.168856, 0.623913, 3.42525, 0.168856, 0.623913, 3.42525]
#K = length(αs)
#rs = [0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 2.53 0 0; 2.53 0 0; 2.53 0 0; -2.53 0 0; -2.53 0 0; -2.53 0 0]
#Rs = [0 0 0; 2.53 0 0; -2.53 0 0]
#Zs = [4 1 1]

S = generate_overlap_matrix(rs, αs)
S_eigen = eigen(S)
X = S_eigen.vectors * diagm(1 ./ sqrt.(S_eigen.values))
H_core = generate_H_core(rs, αs, Rs, Zs)
P_0 = eigen(H_core).vectors

two_electron_ints = generate_two_electron_integrals(rs, αs)
generate_G_STO_nG(P) = generate_G_from_ints(P, two_electron_ints)

C, orbital_energies, electronic_energy = SCF_iteration(N, P_0, H_core, generate_G_STO_nG, X, 1e-6)
total_energy = electronic_energy + nuclear_energy(Rs, Zs)

orbitals = get_orbitals(C, αs, rs)
plot_N = 400
plot_xs = Array(range(-5, 5, plot_N))
plot_ys = Array(range(-5, 5, plot_N))
plot_zs = zeros(plot_N)
println()
println("electronic energy: ", electronic_energy)
println("total energy (with nucleus-nucleus interaction): ", total_energy)
wavefunction = orbitals[3].(plot_xs', plot_ys, plot_zs')
gr()
probrange = maximum(abs.(wavefunction).^2)
orb = heatmap(plot_xs, plot_ys, abs.(wavefunction).^2, aspect_ratio=:equal, grid=false, axis=false, color=cgrad(:inferno, scale=:exp), rightmargin=10Plots.mm, leftmargin=-30Plots.mm, clims=(0, probrange), dpi=1000,xtickfontsize=12,ytickfontsize=12,xguidefontsize=12,yguidefontsize=12,legendfontsize=12)
waverange = maximum(abs.(wavefunction))