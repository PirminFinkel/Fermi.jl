using Plots

"""-------------------------------------
Plot the functions 
-------------------------------------"""

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
        worstcase[i] = min(prod(localdims[1:i]), prod(localdims[i:end]))
    end
    return worstcase
end
    
function plot_ranks(output::String, ranks::Vector{Int64}, I::Int64)
    tick = 13
    guide = 13

    worstcase = worstcase(length(ranks),I)
    p = plot(1:length(localdims)-1, [worstcase, ranks], label = ["worstcase" "TCI ranks"],  yscale= :log10, legendfontsize=guide, xtickfontsize=tick, ytickfontsize=tick, xguidefontsize=guide, yguidefontsize=guide)
    title!("TT rank of A, with R = $R, tol = $tol")
    xlabel!("TT chain position")
    ylabel!("rank")
    savefig(p, output) 
end

function moleculePlot()
    Hpoints = splitH(input)[2]
    Cpoints = splitH(input)[1]
    Opoints = splitH(input)[3]
#SBATCH --time=24:0:0

    plt3d = scatter(Opoints[:,1], Opoints[:,2], Opoints[:,3], seriestype=:scatter, markersize = 7, color = "red", label = "Oxygen")
    scatter!(Hpoints[:,1],Hpoints[:,2], Hpoints[:,3], markersize = 3, color = "blue", label = "Hydrogen")
    scatter!(Cpoints[:,1],Cpoints[:,2], Cpoints[:,3], markersize = 5, color = "yellow", label = "Carbon")
    xlims!(xmin[1], xmax[1])
    ylims!(xmin[2], xmax[2])
    zlims!(xmin[3], xmax[3])

    display(plt3d)
end

function atomPlot(BS::BasisSet, I::Int64)
    limit = 1e-16
    i,li,A,R1,l1 = parametersPhi(BS, I)
    #R1 = round.(R1, digits=2)
    xaxis1 = range(-3, 3, length = 1000)
    xaxis2 = range(xmin[1], xmax[1], length = 1000)
    ydata1 = [abs(phi(A,l1,[k,0.,0.]))>limit ? abs(phi(A,l1,[k,0.,0.])) : limit for k in xaxis1]
    ydata2 = [abs(phi(A,l1,[k-R1[1],0.,0.]))>limit ? abs(phi(A,l1,[k-R1[1],0.,0.])) : limit for k in xaxis2]
    ydata3 = [phi(A,l1,[k,0.,0.]) for k in xaxis1]
    ydata4 = [phi(A,l1,[k-R1[1],0.,0.]) for k in xaxis2]
    plt1 = plot(xaxis1, ydata3, label = "Basis function $i", xlabel = "r in angstrom", ylabel = "φ(r)")
    plt2 = plot(xaxis1, abs.(ydata1), label = "Basis function $i",yscale = :log10, xlabel = "r in angstrom", ylabel = "φ(r)")
    plt3 = plot(xaxis2, ydata4, label = "Basis function $i", xlabel = "r in angstrom", ylabel = "φ(r)")
    plt4 = plot(xaxis2, abs.(ydata2), label = "Basis function $i",yscale = :log10, xlabel = "r in angstrom", ylabel = "φ(r)")
    plt = plot(plt1,plt2,plt3,plt4, layout=(2,2), legend=false, dpi = 600)
    display(plt)
    #plt = plot(plt1,plt2, layout=(1,2), legend=false, dpi = 600, title = "Basisfunction at $R1 with l_xyz = $l1",  titlefont = font(12))
end

function twoAtomPlot(BS::BasisSet, I::Int64, J::Int64)
    i,li,A,R1,l1 = parametersPhi(BS, I)
    j,lj,B,R2,l2 = parametersPhi(BS, J)
    #R1 = round.(R1, digits=2)
    step = nor#SBATCH --time=24:0:0
    [phi(A,l1, (R2-R1)*k/1000) for k in -1000:1:2000]
    ydata3 = [phi(B,l2,(R1-R2+(R2-R1)*k/1000)) for k in -1000:1:2000]
    plt1 = plot(xaxis1,[ydata2,ydata3], label = ["Basis function $i","Basis function $j"], xlabel = "r in angstrom", ylabel = "φ(r)")
    plt2 = plot(xaxis1,[norm.(ydata2),norm.(ydata3)], label = ["Basis function $i","Basis function $j"],yscale =:log10, xlabel = "r in angstrom", ylabel = "φ(r)")
    plt3 = plot(xaxis1, ydata1, label = "overlap", xlabel = "r in angstrom", ylabel = "φ(r)")
    plt4 = plot(xaxis1, norm.(ydata1), label = "overlap",yscale = :log10, xlabel = "r in angstrom", ylabel = "φ(r)")
    plt = plot(plt1,plt3, layout=(1,2), legend=false, dpi = 600)
    display(plt)
    #plt = plot(plt1,plt2, layout=(1,2), legend=false, dpi = 600, title = "Basisfunction at $R1 with l_xyz = $l1",  titlefont = font(12))
end

function generate_Plots(json_filename::String)
    dict = read_json(json_filename)
    for (key,value) in dict
        BS = BasisSetConstructor(key[1:end-2]*".txt", parse(Int64, key[end-1:end]))
        println(key, BS.nbas)
        #plot_ranks(key, value[2],BS.nbas )
    end
end