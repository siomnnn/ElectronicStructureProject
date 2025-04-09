using Plots
using LaTeXStrings

# (values obtained by copy-pasting terminal output from numericalhydrogenplots.jl... then hand-tweaked plots to values)

# for n=1, l=0, m=0
Es_100 = [-0.3527459287219774, -0.44144715216734176, -0.4823846711558115, -0.49533835508703505, -0.49881586073561296, -0.4997026216531426, -0.49992555274058037]
abs_errs_100 = [0.012108411504866232, 0.012067381767884439, 0.010455847929653983, 0.004439605907532396, 0.0014201175137853461, 0.0003996688098742984, 0.00010585944610908005]
rel_errs_100 = [184.02493916487595, 6.380880291848204, 0.8188110953983593, 0.16943760984015685, 0.03767662408424709, 0.010618003517388942, 0.23379717500271005]
l2_errs_100 = [0.10221438472405976, 0.00925714464869682, 0.0005700844962656342, 3.371638536458996e-5, 5.629577240456558e-6, 3.97589433758108e-7, 9.86343989534772e-9]
E_errs_100 = abs.(Es_100.+0.5)

# for n=2, l=1, m=0
Es_210 = [-0.13300624261803295, -0.12671010447849773, -0.1254034731816807, -0.12509953017734166, -0.12502482290495576, -0.12500620501533105, -0.12500155158776738]
abs_errs_210 = [0.020403560759254216, 0.005521251375349444, 0.0012221679679828498, 0.0002962647650837716, 7.349695510873444e-5, 1.8387811002149035e-5, 4.5959894216865416e-6]
rel_errs_210 = [0.5080870442807811, 0.13537509129676442, 0.03361095608337791, 0.010399509685787534, 0.0031265169227789815, 0.0009167466714711905, 0.0003293390179258298]
l2_errs_210 = [0.007353547919257356, 0.0003669286301882497, 2.0438663420353185e-5, 3.3531968141303764e-7, 1.734897487514847e-8, 1.0452429520022224e-9, 9.040026058163646e-11]
E_errs_210 = abs.(Es_210.+0.125)

# for n=3, l=1, m=0
Es_310 = [-0.0592987951344331, -0.05633785902300836, -0.055741046572211565, -0.055601398831280235, -0.05556699462641479, -0.05555841535852779, -0.05555627068094805]
abs_errs_310 = [0.014114732381080491, 0.0036094780759842043, 0.0007966502044322865, 0.0001932161201513051, 4.821069048266402e-5, 1.2030535598214764e-5, 3.00275328755234e-6]
rel_errs_310 = [0.68039934190823, 0.2798802879450375, 0.18699047987143116, 0.10720307172964506, 0.09804982093389725, 0.0957555500785519, 0.0949365512074652]
l2_errs_310 = [0.022315784345665956, 0.0010339150071459655, 5.84060809672424e-5, 3.8089874371131525e-6, 1.5235014438310282e-7, 9.528081947772695e-9, 5.845168969480378e-10]
E_errs_310 = abs.(Es_310.+1/18)

# for n=4, l=3, m=1
Es_431 = [-0.031283121718706355, -0.03125838959576469, -0.03125211258855657, -0.03125053013569312, -0.031250132788051666, -0.031250033229818175, -0.03125000831768803]
abs_errs_431 = [0.0006953895318284668, 0.00034978006468212674, 0.00017579697945608973, 8.807071285448283e-5, 1.714330877162053e-5, 7.3245178432089736e-6, 2.164850423404252e-8]
rel_errs_431 = [5.1819031841441054, 6.34835357591059, 6.623825953961949, 6.332021943433294, 0.930873456772088, 0.9255825600711227, 0.9216509213122697]
l2_errs_431 = [2.3349412224116366e-5, 9.242848979959859e-5, 1.6504998897452177e-5, 1.2375605071653902e-6, 2.725949393708823e-6, 4.988751470805399e-7, 5.8975587375722445e-12]
E_errs_431 = abs.(Es_431.+0.03125)

pyplot()
plot(Ns, 256*Ns.^(-2), xscale=:log2, yscale=:log2, xlims=(2^6, 2^12), xlabel=L"$N$ [-]", ylabel="error [-]", legend=:bottomleft, color=:black, linestyle=:dash, primary=false)
plot!(Ns, 1/256*Ns.^(-2), color=:black, linestyle=:dash, primary=false)
plot!(Ns, Ns.^(-2), color=:black, linestyle=:dash, primary=false)
plot!(Ns, abs_errs_310, label="max. absolute error", color=:firebrick)
plot!(Ns, l2_errs_310, label=L"$L^2$ error", color=:royalblue)
plot!(Ns, E_errs_310, label="energy error", color=:orange)
plot!([], [], label=L"\mathcal{O}(h^2)", color=:black, linestyle=:dot)
savefig("plots/numerical/errors310.png")

plot(Ns, 4096*Ns.^(-2), xscale=:log2, yscale=:log2, xlims=(2^6, 2^12), xlabel=L"$N$ [-]", ylabel="error [-]", legend=:bottomleft, color=:black, linestyle=:dash, primary=false)
plot!(Ns, 1/16*Ns.^(-2), color=:black, linestyle=:dash, primary=false)
plot!(Ns, 16*Ns.^(-2), color=:black, linestyle=:dash, primary=false)
plot!(Ns, abs_errs_100, label="max. absolute error", color=:firebrick)
plot!(Ns, l2_errs_100, label=L"$L^2$ error", color=:royalblue)
plot!(Ns, E_errs_100, label="energy error", color=:orange)
plot!([], [], label=L"\mathcal{O}(h^2)", color=:black, linestyle=:dot)
savefig("plots/numerical/errors100.png")

plot(Ns, Ns.^(-2), xscale=:log2, yscale=:log2, xlims=(2^6, 2^12), xlabel=L"$N$ [-]", ylabel="error [-]", legend=:bottomleft, color=:black, linestyle=:dash, primary=false)
plot!(Ns, 1/256*Ns.^(-2), color=:black, linestyle=:dash, primary=false)
plot!(Ns, 256Ns.^(-2), color=:black, linestyle=:dash, primary=false)
plot!(Ns, abs_errs_431, label="max. absolute error", color=:firebrick)
plot!(Ns, l2_errs_431, label=L"$L^2$ error", color=:royalblue)
plot!(Ns, E_errs_431, label="energy error", color=:orange)
plot!([], [], label=L"\mathcal{O}(h^2)", color=:black, linestyle=:dot)
savefig("plots/numerical/errors431.png")

plot(Ns, 4096*Ns.^(-2), xscale=:log2, yscale=:log2, xlims=(2^6, 2^12), xlabel=L"$N$ [-]", ylabel="error [-]", legend=:topright, color=:black, linestyle=:dash, primary=false)
plot!(Ns, 1/1024*Ns.^(-2), color=:black, linestyle=:dash, primary=false)
plot!(Ns, 4*Ns.^(-2), color=:black, linestyle=:dash, primary=false)
plot!(Ns, 1024*4096*Ns.^(-2), color=:black, linestyle=:dash, primary=false)
plot!(Ns, abs_errs_430, label="max. absolute error", color=:firebrick)
plot!(Ns, l2_errs_430, label=L"$L^2$ error", color=:royalblue)
plot!(Ns, E_errs_430, label="energy error", color=:orange)
plot!([], [], label=L"\mathcal{O}(h^2)", color=:black, linestyle=:dot)
savefig("plots/numerical/errors430.png")

plot(Ns, 4096*Ns.^(-2), xscale=:log2, yscale=:log2, xlims=(2^6, 2^12), xlabel=L"$N$ [-]", ylabel="error [-]", legend=:topright, color=:black, linestyle=:dash, primary=false)
plot!(Ns, 1/16*Ns.^(-2), color=:black, linestyle=:dash, primary=false)
plot!(Ns, 16*Ns.^(-2), color=:black, linestyle=:dash, primary=false)
plot!(Ns[3:end], 256*4096*Ns[3:end].^(-2), color=:black, linestyle=:dash, primary=false)
plot!(Ns, old_abs_errs_100, label="max. absolute error", color=:firebrick)
plot!(Ns, old_l2_errs_100, label=L"$L^2$ error", color=:royalblue)
plot!(Ns, old_E_errs_100, label="energy error", color=:orange)
plot!([], [], label=L"\mathcal{O}(h^2)", color=:black, linestyle=:dot)
savefig("plots/numerical/old_errors100.png")


