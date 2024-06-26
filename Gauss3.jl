using Fermi
using GaussianBasis 
using StaticArrays
using Plots
using LinearAlgebra
using QuanticsTCI
import TensorCrossInterpolation as TCI
using Profile
using StatProfilerHTML
using TCIITensorConversion
using ITensors

include("setup.jl")

const ang2bohr = 1.8897261246257702

"""-------------------------
Here we define our functions
-------------------------"""
#Basis function
function phi(BS::BasisSet, numberBasisfunction::Int64, l_xyz::Int64, r::Vector{Float64})
    R = BS.basis[numberBasisfunction].atom.xyz
    r = (r-R).*ang2bohr

    basisfunction = BS.basis[numberBasisfunction]
    L = permute(basisfunction.l)
    l = L[l_xyz,:]

    result = 0
    for i in 1:length(basisfunction.coef)
        result += basisfunction.coef[i]*exp(-basisfunction.exp[i]*(r[1]^2+r[2]^2+r[3]^2))
    end
    
    result *= (r[1]^l[1])*(r[2]^l[2])*(r[3]^l[3])
    return result 
end


function phi(coef::AbstractArray, l_xyz::Vector{Int64}, r::AbstractVector)
    r = r.*ang2bohr
    result = 0.
    for i in 1:(size(coef,1))
        result += coef[i,1]*exp(-coef[i,2]*(r[1]^2+r[2]^2+r[3]^2))
    end
    result *= (r[1]^l_xyz[1])*(r[2]^l_xyz[2])*(r[3]^l_xyz[3])
    return result 
end

#Analytic integral
function a(BS::BasisSet, r::Vector{Float64}, i::Int64, j::Int64)
    fakeatom = [GaussianBasis.Atom(-1, 1, r)]
    result = Coulomb_Integral(BS, i, j, fakeatom)
    return result
end

function Coulomb_Integral(BS::BasisSet, i::Int64, j::Int64, atom::Vector{<:Atom}) #choose two Basis functions from the set 
    out = zeros(GaussianBasis.num_basis(BS.basis[i]), GaussianBasis.num_basis(BS.basis[j]))
    GaussianBasis.generate_V_pair!(out, BS.basis[i], BS.basis[j], atom)
    return out
end

#Full Discretized analytic Tensor
function A(c::Config, list::Vector{Int64})
    I = list[1]
    coord = QG.quantics_to_origcoord(c.grid, list[2:end])

    i, li = getAtomL(c.BS, splitIndex(c.BS, I)[1])
    j, lj = getAtomL(c.BS, splitIndex(c.BS, I)[2])

    return a(c.BS, [coord[1], coord[2], coord[3]], i, j)[li, lj]
end

#Double phi
function P(c::Config, list::Vector{Int64})
    I = list[1]

    i,li,A,R1,l1 = parametersPhi(c.BS, splitIndex(c.BS, I)[1])
    j,lj,B,R2,l2 = parametersPhi(c.BS, splitIndex(c.BS, I)[2])

    return phi(A,l1,collect(QG.quantics_to_origcoord(c.grid, list[2:end]))-R1)*phi(B,l2,collect(QG.quantics_to_origcoord(c.grid, list[2:end]))-R2)
end

#Total Coulomb Integral (μν|λρ)
function totalIntegral(c::Config, I,J)

    R = c.grid.R
    i,li,A1,R1,l1 = parametersPhi(c.BS, splitIndex(c.BS, I)[1])
    j,lj,B,R2,l2 = parametersPhi(c.BS, splitIndex(c.BS, I)[2])
    
    result = 0.
    for coord in Iterators.product(1:8, ntuple(x -> 1:8, R - 1)...)
        x = collect(QG.quantics_to_origcoord(c.grid,collect(coord)))
        y = [J,coord...]
        result += phi(A1,l1,x-R1)*phi(B,l2,x-R2)*A(c,y)
    end
    
    return result*prod(collect(QG.grid_max(c.grid))-collect(QG.grid_min(c.grid)))*ang2bohr^3/(2^(3*R))
end

#General Kernel and Basis functions
function generalBasis(grid::QG.Grid, list::Vector{Int64})::Float64
    I,J = splitIndex(c.BS, list[1])
    x = collect(QG.quantics_to_origcoord(grid,collect(list[2:end])))
    return phiGeneral(I,x)*phiGeneral(J,x)
end

function phiGeneral(I::Int64, x::Vector)::Float64
    #to be further defined
    return phi(A,l1,x-R1)
end

function generalKernel(list::Vector{Int64})::Float64
    x = collect(QG.quantics_to_origcoord(grid,collect(list[1:R])))
    y = collect(QG.quantics_to_origcoord(grid2,collect(list[R+1:2*R])))
    return u(x,y)
end

function u(x::Vector{Float64}, y::Vector{Float64})
    #return (x != y ? 1/norm((x-y).*ang2bohr) : 1)
    return 1/norm(x-y).*ang2bohr
end



"""--------------------------------------
print out information
--------------------------------------"""

function infoGeneral(c::Config)
    input = c.name*".txt"
    println("Inputfile: \t"*c.name*".txt","\n")

    println("Number of atoms:\t\t", read_xyz(c.name*".txt")[2])
    println("Number of basis shells:\t\t", length(c.BS.basis))
    println("Number of basis functions: \t", lengthIndex(c.BS))
    println("Length of multiindex I: \t", lengthIndex(c.BS)^2,"\n")

    xmin_molecule = getMaxima(toArray(input))[1]
    xmax_molecule = getMaxima(toArray(input))[2]
    println("Molecule size:\n",Tuple(xmax_molecule) , "\n", Tuple(xmin_molecule), "\n")
    println("Grid size:\n", QG.grid_max(c.grid), "\n", QG.grid_min(c.grid), "\n")
end

function infoBasis(c::Config)
    println("-------------------------------------------------------------------------------------")
    println("I\ti\tAtom\tx\t\ty\t\tz\t\tl\tl_xyz")
    println("-------------------------------------------------------------------------------------")
    currentAtom = Vector{Float64}
    for k in 1:lengthIndex(c.BS)
        i,li = getAtomL(c.BS, k)
        bsf = c.BS.basis[i]
        lxyz = permute(bsf.l)[li,:]
        
        print(k,"\t",i , "\t")
        if bsf.atom.xyz != currentAtom
            temp = round.(bsf.atom.xyz, digits = 4)
            print( GaussianBasis.symbol(bsf.atom), "\t", temp[1],"\t\t", temp[2],"\t\t", temp[3],"\t\t")
            currentAtom = bsf.atom.xyz
        else
            print( "\t"^7)
        end
        println(bsf.l,"\t$(lxyz[1])$(lxyz[2])$(lxyz[3])")

    end
    println("-------------------------------------------------------------------------------------")
end

"""--------------------------------------------------
Tensor cross interpolation to compress the functions 
--------------------------------------------------"""

function worstcase1(R::Int64)
    worstcase = Vector{Float64}(undef,R)
    for i in 1:Int(ceil(R/2))
        worstcase[i] = 8^i
        worstcase[R+1-i] = 8^i
    end
    worstcase
end

function worstcase2(R::Int64, lengthFirst::Int64)
    localdims = fill(8, R+1)
    localdims[1] = lengthFirst^2

    worstcase = Vector{Float64}(undef,R)
    for i in 1:R
        worstcase[i] = min(prod(localdims[1:i]), prod(localdims[i:end]))
    end
    return worstcase
end

function fullCompress(c::Config, tol::Float64)
    R = c.grid.R
    localdims = fill(8, R+1)
    localdims[1] = lengthIndex(c.BS)^2

    A_comp(list::Vector) = A(c, list)

    cA = TensorCrossInterpolation.CachedFunction{Float64}(A_comp, localdims)
    tci, ranks, errors = crossinterpolate2(Float64, cA, localdims; tolerance = tol)
    return r = TensorCrossInterpolation.linkdims(tci)

    worstcase = worstcase2(R, c.BS.nbas)

    tick = 13
    guide = 13

    p = plot(1:length(localdims)-1, [worstcase, r], label = ["worstcase" "TCI ranks"],  yscale= :log10, legendfontsize=guide, xtickfontsize=tick, ytickfontsize=tick, xguidefontsize=guide, yguidefontsize=guide)
    title!("TT rank of A, with R = $R, tol = $tol")
    xlabel!("TT chain position")
    ylabel!("rank")
    savefig(p, generateOutputName(c.name, R, tol))
end

function makeTT(BS::BasisSet)
    localdims = fill(8, R+1)
    localdims[1] = lengthIndex(BS)^2
    
    #Tensor A
    cA = TensorCrossInterpolation.CachedFunction{Float64}(A, localdims)
    tci_A, ranks_A, errors_A = crossinterpolate2(Float64, cA, localdims; tolerance = tol)

    #Tensor P 
    cP = TensorCrossInterpolation.CachedFunction{Float64}(P, localdims)
    tci_P, ranks_P, errors_P = crossinterpolate2(Float64, cP, localdims; tolerance = tol)


    return (tci_A, tci_P)
end

function makeGeneralTT(A_general::Function, B_general::Function, C_general::Function, basislength::Int64)

    localdims = fill(8, R+1)
    localdims[1] = basislength^2
    kerneldims = fill(8,2*R)

    #Tensor A
    cA = TensorCrossInterpolation.CachedFunction{Float64}(A_general, localdims)
    tci_A, ranks_A, errors_A = crossinterpolate2(Float64, cA, localdims; tolerance = tol)

    #Tensor B 
    cB = TensorCrossInterpolation.CachedFunction{Float64}(B_general, kerneldims)
    tci_B, ranks_B, errors_B = crossinterpolate2(Float64, cB, kerneldims; tolerance = tol)

    #Tensor C
    cC = TensorCrossInterpolation.CachedFunction{Float64}(C_general, localdims)
    tci_C, ranks_C, errors_C = crossinterpolate2(Float64, cC, localdims; tolerance = tol)

    return (tci_A, tci_B, tci_C)
end

function evaluateTT(tci_A, tci_P)
    tt1 = TCI.TensorTrain(tci_A)
    M1 = ITensors.MPS(tt1)

    tt2 = TCI.TensorTrain(tci_P)
    M2 = ITensors.MPS(tt2)

    for k in 2:R+1
        replaceind!(M1[k],inds(M1[k])[2], inds(M2[k])[2])
    end

    #M = ITensors.contract(M1)*ITensors.contract(M2)

    result = M1[R+1]*M2[R+1]
    for k in R:-1:1
        result = M1[k]*(M2[k]*result)
    end

    return Matrix(result, inds(result)[1], inds(result)[2])*norm(prod(xmax-xmin)*ang2bohr^3)/(2^(3*R))
    #return Matrix(M, inds(M)[1], inds(M)[2])*norm(prod(xmax-xmin)*ang2bohr^3)/(2^(3*R))
end

function evaluateGeneralTT(tci_A, tci_B, tci_C)
    mpsA = ITensors.MPS(TCI.TensorTrain(tci_A))
    mpsB = ITensors.MPS(TCI.TensorTrain(tci_B))
    mpsC= ITensors.MPS(TCI.TensorTrain(tci_C))

    for k in 2:R+1
        temp = (k==2  ? 1 : 2)
        replaceind!(mpsB[k-1],inds(mpsB[k-1])[temp], inds(mpsA[k])[2])
        replaceind!(mpsB[R+k-1],inds(mpsB[R+k-1])[2], inds(mpsC[k])[2])
    end

    #return mpsA, mpsB, mpsC

    #lower contraction
    result = mpsC[R+1]*mpsB[2*R]
    for k in R:-1:2
        result = mpsC[k]*(mpsB[R+k-1]*result)
    end
    result = mpsC[1]*result

    #upper contraction
    result = mpsA[R+1]*(mpsB[R]*result)
    for k in R:-1:2
        result = mpsA[k]*(mpsB[k-1]*result)
    end
    result = mpsA[1]*result

    #M = ITensors.contract(ITA)*ITensors.contract(ITB)*ITensors.contract(ITC)
    return Matrix(result, inds(result)[1], inds(result)[2])*norm(prod(xmax-xmin)*ang2bohr^3)/(2^(3*R))
end



"""-------------------------------------
Plot the functions 
-------------------------------------"""

function moleculePlot()
    Hpoints = splitH(input)[2]
    Cpoints = splitH(input)[1]
    Opoints = splitH(input)[3]

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
    step = norm((R1-R2)/1000)*ang2bohr
    xaxis1 = range(-1000*step, 2000*step, length = 3001)
    ydata1 = [phi(A,l1,(R2-R1)*k/1000)*phi(B,l2,(R1-R2+(R2-R1)*k/1000)) for k in -1000:1:2000]
    ydata2 = [phi(A,l1, (R2-R1)*k/1000) for k in -1000:1:2000]
    ydata3 = [phi(B,l2,(R1-R2+(R2-R1)*k/1000)) for k in -1000:1:2000]
    plt1 = plot(xaxis1,[ydata2,ydata3], label = ["Basis function $i","Basis function $j"], xlabel = "r in angstrom", ylabel = "φ(r)")
    plt2 = plot(xaxis1,[norm.(ydata2),norm.(ydata3)], label = ["Basis function $i","Basis function $j"],yscale =:log10, xlabel = "r in angstrom", ylabel = "φ(r)")
    plt3 = plot(xaxis1, ydata1, label = "overlap", xlabel = "r in angstrom", ylabel = "φ(r)")
    plt4 = plot(xaxis1, norm.(ydata1), label = "overlap",yscale = :log10, xlabel = "r in angstrom", ylabel = "φ(r)")
    plt = plot(plt1,plt3, layout=(1,2), legend=false, dpi = 600)
    display(plt)
    #plt = plot(plt1,plt2, layout=(1,2), legend=false, dpi = 600, title = "Basisfunction at $R1 with l_xyz = $l1",  titlefont = font(12))
end

"""-------------------------------------------------------------
We define a test function, to check our functions on correctness
-------------------------------------------------------------"""

function test_A(c::Config)
    R = c.grid.R
    I = 1 
    N = 20
    #A_tci = fullCompress()[1]

    i,li,A,R1,l1 = parametersPhi(c.BS, splitIndex(c.BS, I)[1])
    j,lj,B,R2,l2 = parametersPhi(c.BS, splitIndex(c.BS, I)[2])

    v = Vector{Float64}(undef, N)
    for k in 1:N

        x = collect(QG.quantics_to_origcoord(c.grid,rand(1:8,R)))
        if k ==1 
            x = collect(QG.quantics_to_origcoord(c.grid,fill(4,R)))
        end

        result = 0.
        for coord in Iterators.product(1:8, ntuple(x -> 1:8, R - 1)...)
            y = collect(QG.quantics_to_origcoord(c.grid,collect(coord)))
            if x != y
                result += phi(A,l1,(y-R1))*phi(B,l2,(y-R2))/(norm(x-y)*ang2bohr)
            end
        end
        
        num = result*prod(collect(QG.grid_max(c.grid))-collect(QG.grid_min(c.grid)))*ang2bohr^3/(2^(3*R))
        ana = a(c.BS,x, i,j)[li,lj]
        
        v[k] = abs.(num-ana)

    end

    compareResults2(sum(v)/N, maximum(v))
    
end

function test_phi(c::Config)
    R = c.grid.R
    I = 1
    N = lengthIndex(c.BS)
    num = zeros(Float64, N,N)
    
    for n in 1:N
        for m in 1:N
            if m>=n 
                i,li,A,R1,l1 = parametersPhi(c.BS, n)
                j,lj,B,R2,l2 = parametersPhi(c.BS, m)

                result = 0.
                for coord in Iterators.product(1:8, ntuple(x -> 1:8, R - 1)...)
                    y = collect(QG.quantics_to_origcoord(c.grid,collect(coord)))
                    result += phi(A,l1,(y-R1))*phi(B,l2,(y-R2))
                end
                num[n,m] = num[m,n] = result*prod(collect(QG.grid_max(c.grid))-collect(QG.grid_min(c.grid)))*ang2bohr^3/(2^(3*R))
            end
        end
    end

    ana = overlap(c.BS, c.BS)
    result = (num-ana)
    compareResults2(sum(result)/prod(size(result)), maximum(result))
end

function compressPhi(c::Config, tol::Float64)
    R = c.grid.R
    I=1
    localdims = fill(8, R)
    N = lengthIndex(c.BS)
    num = zeros(Float64, N,N)
    
    for n in 1:N
        for m in 1:N
            if m>=n 

                i,li,A,R1,l1 = parametersPhi(c.BS, n)
                j,lj,B,R2,l2 = parametersPhi(c.BS, m)

                p(list::Vector) = phi(A,l1,collect(QG.quantics_to_origcoord(c.grid, list))-R1)*phi(B,l2,collect(QG.quantics_to_origcoord(c.grid, list))-R2)

                cp = TensorCrossInterpolation.CachedFunction{Float64}(p, localdims)
                pivot3 = collect(QG. origcoord_to_quantics(c.grid,Tuple(Float64(x) for x in (R1+R2)/2)))
                tci = crossinterpolate2(Float64, cp, localdims, [pivot3]; tolerance = tol)[1]

                num[n,m] = num[m,n] = sum(tci)*prod(collect(QG.grid_max(c.grid))-collect(QG.grid_min(c.grid)))*ang2bohr^3/(2^(3*R))
            end
        end
    end

    ana = overlap(c.BS, c.BS)
    return num-ana
end

function test_total(c::Config)
    N = 20
    BS = c.BS
    L = lengthIndex(BS)^2

    v = Vector{Float64}(undef, N)     
    for k in 1:N
        I,J = rand(1:L,2)

        i, li = getAtomL(BS, splitIndex(BS, I)[1])
        j, lj = getAtomL(BS, splitIndex(BS, I)[2])
        n, ln = getAtomL(BS, splitIndex(BS, J)[1])
        m, lm = getAtomL(BS, splitIndex(BS, J)[2])

        v[k] =  abs( totalIntegral(c,I,J)-ERI_2e4c(BS,i,j,n,m)[li,lj,ln,lm])
    end
    compareResults2(sum(v)/N, maximum(v))
end

function test_P(c::Config)
    R = c.grid.R
    I = 6
    N = 20
    #A_tci = fullCompress()[1]

    i,li,A,R1,l1 = parametersPhi(c.BS, splitIndex(c.BS, I)[1])
    j,lj,B,R2,l2 = parametersPhi(c.BS, splitIndex(c.BS, I)[2])

    v = Vector{Float64}(undef,N)
    for k in 1:N
        y = rand(1:8,R)
        x = collect(QG.quantics_to_origcoord(c.grid,y))

        result = phi(A,l1,(x-R1))*phi(B,l2,(x-R2))
        ana = P(c,[I,y...])
        
        v[k] = abs(result-ana)
    end
    compareResults2(sum(v)/N, maximum(v))
end

function generateAna(BS::BasisSet)
    N = lengthIndex(BS)^2
    Arr = Matrix{Float64}(undef,N,N)
    C = ERI_2e4c(BS)

    for k1 in 1:N
        for k2 in 1:N
            
        i, j = splitIndex(c.BS, k1)
        n, m = splitIndex(c.BS, k2)

        Arr[k1,k2] = C[i,j,n,m]
        end
    end
    return Arr
end

"-----------------------------------------------------------------------------------------------------------------------------"

function main()
    mol1 = ("hydrogen","water","ethane", "propanol", "glucose")
    mol = ("hydrogen",)
    v = Vector{Config}(undef,4*length(mol))
    for n in eachindex(mol)
        for m in (1,2)
            v[4*(n-1)+m] =  Config(mol[n]*".txt",m, 5)
            v[4*(n-1)+m+2] = Config(mol[n]*".txt",m,7)
        end
    end

    for n in eachindex(v)
        c = v[n]
        t= @elapsed ranks = fullCompress(c,1e-7)
        println(ranks ,t)
    end
end