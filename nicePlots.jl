using Plots

"""--------------------------------------------------------------------------------
Ranks and JSON 
--------------------------------------------------------------------------------"""

function worstcase(R::Int64)
    worstcase = Vector{Float64}(undef,R)
    for i in 1:Int(ceil(R/2))
        worstcase[i] = 8^i
        worstcase[R+1-i] = 8^i
    end
    worstcase
end

function worstcase(R::Int64, I::Int64)
    localdims = fill(8, R+1)
    localdims[1] = I^2

    worstcase = Vector{Float64}(undef,R)
    for i in 1:R
        worstcase[i] = min(prod(localdims[1:i]), prod(localdims[i+1:end]))
    end
    return worstcase
end

function plot_ranks(output::String, ranks::Vector{Int64}, I::Int64)
    tick = 13
    guide = 13

    worst = worstcase(length(ranks),I)
    p = plot(1:length(ranks), [worst, ranks], label = ["worstcase" "TCI ranks"],  yscale= :log10, legendfontsize=guide, xtickfontsize=tick, ytickfontsize=tick, xguidefontsize=guide, yguidefontsize=guide)
    title!("TT rank of A, with R = 5, tol = 1e-7")
    xlabel!("TT chain position")
    ylabel!("rank")
    display(p)
    savefig(p, "./plots/"*output) 
end

function plot_ranks(output::String, ranks::Vector{Vector{Int64}}, I::Int64)
    tick = 13
    guide = 13

    worst = worstcase(length(ranks[1]),I)
    p = plot(1:length(ranks[1]), [worst, ranks...], label = ["worstcase" "tol = 1e-7" "tol = 1e-5" "tol = 1e-3" ],  yscale= :log10, legendfontsize=guide, xtickfontsize=tick, ytickfontsize=tick, xguidefontsize=guide, yguidefontsize=guide)
    title!("TT rank of A, with R = 5")
    xlabel!("TT chain position")
    ylabel!("rank")
    display(p)
    savefig(p, "./plots/"*output) 
end

function plot_json(json_path::String)
    dict = read_json(json_path)

    for (name, value) in dict
        plot_ranks(name*"-3", Int.(value[2]), value[3])
    end

end

function plot_combined_json(json_path1::String, json_path2::String, json_path3::String)
    dict1 = read_json(json_path1)
    dict2 = read_json(json_path2)
    dict3 = read_json(json_path3)

    for (name, value) in dict3
        plot_ranks(name*"_combined", [Int.(value[2]), Int.(dict1[name][2]), Int.(dict2[name][2])], value[3])
    end
end


"""--------------------------------------------------------------------------------
New Plots
--------------------------------------------------------------------------------"""

function plot_P(I::Int64, J::Int64, BS::BasisSet)
    temp1, temp2, A,R1,l1 = parametersPhi(BS, I)
    temp3, temp4, B,R2,l2 = parametersPhi(BS, J)    

    guide = tick = 11
    s = (800,800)
    
    dir = Vector{Float64}(undef,3)
    xdata = range(-1.5, 1.5, length = 300)

    if R1!=R2 
        dir = (R2-R1)/norm(R2-R1)
        dist = norm(R2-R1)
        xdata = range(-dist, 2*dist, length = 300)
    else 
        dir = [1.,0.,0.]
    end

    overlap = [phi(A,l1,k*dir)*phi(B,l2,R2-(R1+k*dir)) for k in xdata]
    y2 = [phi(B,l2,R2-(R1+k*dir)) for k in xdata]
    y1 = [phi(A,l1,R1-(R1+k*dir)) for k in xdata]

    plot(xdata, [y1 y2] ,label = [ "phi1" "phi2"], xlabel= "r-R1 in angstrom", ylabel="P(r)", legend = true, legendfontsize=guide, xtickfontsize=tick, ytickfontsize=tick, xguidefontsize=guide, yguidefontsize=guide, size=s)
    plot!(xdata, overlap ,label = "overlap",linewidth = 2, color = :black)
    return current()
end

function plot_A(I::Int64, J::Int64, BS::BasisSet)
    guide = tick = 11
    s = (800,800)

    i, li = getAtomL(BS, I)
    j, lj = getAtomL(BS, J)
    R1 = BS.basis[i].atom.xyz
    R2 = BS.basis[j].atom.xyz

    f(r::Vector{Float64}) = a(BS, r, i, j)[li, lj]

    dir = [1.,0.,0.]
    xdata = range(-10., 20., length = 300)

    if R1!=R2 
        dir = (R2-R1)/norm(R2-R1)
        dist = norm(R2-R1)*10
        xdata = range(-dist, 2*dist, length = 300)
    end

    overlap = [ a(BS, R1+k*dir, i, j)[li, lj] for k in xdata]

    plot(xlabel= "r-R1 in angstrom", ylabel="P(r)", legendfontsize=guide, tickfontsize=tick,  guidefontsize=guide, size=s)
    plot!(xdata, overlap ,label = "A(r)")

    temp1, temp2, A,R1,l1 = parametersPhi(BS, I)
    temp3, temp4, B,R2,l2 = parametersPhi(BS, J)    
    overlap2 = [phi(A,l1,k*dir)*phi(B,l2,R2-(R1+k*dir)) for k in xdata] #/abs((2-k)*ang2bohr)
    plot!(xdata, overlap2 ,label = "overlap",linewidth = 2, color = :black)

    hyperbel = [1/abs(20-k) for k in xdata]
    #plot!(xdata, hyperbel ,label = "coulomb kernel")
    return current()
end


function plot_function_plane(f::Function, Offset)
    tick = guide = 11
    Res = 300
    s = (1500, 450)
    incr = 0.1
    tol = 10e-7
     
    xs = range(Offset[1,1], Offset[1,2], Res)
    ys = range(Offset[2,1], Offset[2,2], Res)
    zs = range(Offset[3,1], Offset[3,2], Res)

    p = Vector{Any}(undef,3)
    p[1] = surface(xs,ys, (x,y) -> f([x,y,0.]),xlabel= "x", ylabel = "y", zlabel="f(x,y,0)", legend = false, legendfontsize=guide, xtickfontsize=tick, ytickfontsize=tick, xguidefontsize=guide, yguidefontsize=guide, size=s)
    p[2] = surface(ys,zs, (y,z) -> f([0.,y,z]),xlabel= "y", ylabel = "z", zlabel="f(0,y,z)", legend = false, legendfontsize=guide, xtickfontsize=tick, ytickfontsize=tick, xguidefontsize=guide, yguidefontsize=guide, size=s)
    p[3] = surface(xs,zs, (x,z) -> f([x,0.,z]),xlabel= "x", ylabel = "z", zlabel="f(x,0,z)", legend = false, legendfontsize=guide, xtickfontsize=tick, ytickfontsize=tick, xguidefontsize=guide, yguidefontsize=guide, size=s)

    #plt = plot(x for x in p , title=["" "Axis Plot for r2 = $rj" ""], layout=(1,3))
    plt = plot(p[1], p[2], p[3], title=["" "Axis plot of $(GaussianBasis.symbol(BS.basis[i].atom))-atom basis function with l = $l1" ""], layout=(1,3))
    return plt
end

function plot_function_plane(BS::BasisSet, I::Int64)
     
    i, li, A,R1,l1 = parametersPhi(BS, I)
    f(r) = phi(A,l1,r)

    Offset = incr
    while (abs(f([Offset, Offset,0.])) > tol) || (abs(f([Offset,Offset, 0.])) > tol) || (abs(f([0.,Offset,Offset])) > tol)
        Offset += incr
    end
       
   return plot_function_plane(f, fill(Offset, (2,3)))
end

function plot_function_axis(BS::BasisSet, I::Int64)
    tick = guide = 11
    Res = 300
    s = (1500, 450)
    incr = 0.1
    tol = 10e-4
    
    i, li, A,R1,l1 = parametersPhi(BS, I)

    f(r) = phi(A,l1,r)

    Offset = incr
    while (abs(f([Offset, 0.,0.])) > tol) || (abs(f([0.,Offset, 0.])) > tol) || (abs(f([0.,0.,Offset])) > tol)
        Offset += incr
    end
    
    x = range(-Offset, +Offset, Res)
    y = range(-Offset, +Offset, Res)
    z = range(-Offset, +Offset, Res)

    if sum(l1) == 0
        s = (800,600)
        plot(x, x -> f([x,0.,0.]), xlabel= "x", ylabel="f(x,0,0)", legend = false, legendfontsize=guide, xtickfontsize=tick, ytickfontsize=tick, xguidefontsize=guide, yguidefontsize=guide, size=s)
        title!("Axis plot of $(GaussianBasis.symbol(BS.basis[i].atom))-atom basis function with l = $l1")
        return current()
    end

    p = Vector{Any}(undef,3)
    p[1] = plot(x, x -> f([x,0.,0.]), xlabel= "x", ylabel="f(x,0,0)", legend = false, legendfontsize=guide, xtickfontsize=tick, ytickfontsize=tick, xguidefontsize=guide, yguidefontsize=guide, size=s)
    p[2] = plot(y, y -> f([0.,y,0.]), xlabel= "y", ylabel="f(0,y,0)", legend = false, legendfontsize=guide, xtickfontsize=tick, ytickfontsize=tick, xguidefontsize=guide, yguidefontsize=guide, size=s)
    p[3] = plot(z, z -> f([0.,0.,z]), xlabel= "z", ylabel="f(0,0,z)", legend = false, legendfontsize=guide, xtickfontsize=tick, ytickfontsize=tick, xguidefontsize=guide, yguidefontsize=guide, size=s)

    #plt = plot(x for x in p , title=["" "Axis Plot for r2 = $rj" ""], layout=(1,3))
    plt = plot(p[1], p[2], p[3], title=["" "Axis plot of $(GaussianBasis.symbol(BS.basis[i].atom))-atom basis function with l = $l1" ""], layout=(1,3))
    return plt
end

function plot_multiple_in_x(BS, list, atom)
    tick = guide = 11
    Res = 300
    s = (800,600)
    incr = 0.1
    tol = 10e-4

    max = 0.0
    for I in list 
        i, li, A,R1,l1 = parametersPhi(BS, I)
        f(r) = phi(A,l1,r)

        Offset = incr
        while (abs(f([Offset, 0.,0.])) > tol) || (abs(f([0.,Offset, 0.])) > tol) || (abs(f([0.,0.,Offset])) > tol)
            Offset += incr
        end

        if Offset > max
            max = Offset
        end
    end
 
    x = range(-max, +max, Res)

    plot(xlabel= "x", ylabel="f(x,0,0)", legendfontsize=guide, tickfontsize=tick, guidefontsize=guide, size=s)
    for I in list 
        i, li, A,R1,l1 = parametersPhi(BS, I)
        f(r) = phi(A,l1,r)

        plot!(x, x -> f([x,0.,0.]), label = "l = $(sum(l1))")
    end

    title!("$atom basis functions")
    return current()
end




