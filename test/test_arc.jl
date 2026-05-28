

# Single-layer solve on the upper semicircle arc dom=[-π/2, π/2].
# RHS g(θ) = cos(2θ); convergence measured via the potential at an off-arc point.
println("\nArc (unit circle, dom=[-π/2,π/2]): single-layer solve, potential at z=2")
let
    γ(s)   = exp(1im * s)
    γs(s)  = 1im * exp(1im * s)
    γss(s) = -exp(1im * s)
    g(θ) = 0.2 * cos(2θ) + 0.5 * sin(3θ)

    Nref  = 512
    Γref  = arc(γ, γs, γss, Nref; dom=[-π/2, π/2])

    solve(Γref) = V(g.(angle.(Γref.γ)), Γref)
    obj(Γ) = sum(Γ.s .* solve(Γ) .* Γ.w)
    aref = obj(Γref)

    for N in Nspan
        Γ   = arc(N; dom=[-π/2, π/2])
        a = obj(Γ)
        err = abs(aref - a)
        printerror(N, err)
    end
end
