using Fermi
using GaussianBasis 
using StaticArrays
using Plots
using LinearAlgebra
using QuanticsTCI
using TensorCrossInterpolation
using Profile
using StatProfilerHTML

include("setup.jl")

"""--------------------------------------------------------------------------------------------
In this section, we define our preset paremeters, that are later used as input for the function
--------------------------------------------------------------------------------------------"""

input = "molecules.txt" 
R = 5
tol = 1e-7
offset = 5

"----------------------------------------------------------------------------------------------"

#Parameters
lib = 1
basisSetGaussian = ("sto-3g", "def2-svp") 
output = generateOutputName(input,0)  
ang2bohr = 1.8897261246257702

#Basis
BS1 = BasisSet(basisSetGaussian[1], read_xyz(input)[lib], spherical=false, lib =:acsint)

#Grid
xmin_molecule = getMaxima(toArray(input))[2]
xmax_molecule = getMaxima(toArray(input))[1]

xmin = xmin_molecule - [offset, offset, offset]
xmax = xmax_molecule + [offset, offset, offset]

xs = range(xmin[1], xmax[1], length = 2^R)
ys = range(xmin[2], xmax[2], length = 2^R)
zs = range(xmin[3], xmax[3], length = 2^R)
grid = QG.DiscretizedGrid{3}(R, (xmin[1], xmin[2], xmin[3]), (xmax[1], xmax[2], xmax[3]))

"""-------------------------------
Print out the important Parameters
-------------------------------"""

function information()
    println("Input file: \t"*input)
    println("Output file: \t"*output)
    println("Basisset:\t"*basisSetGaussian[lib]*"\n")
   
    println("Number of atoms:\t\t", read_xyz(input)[2])
    println("Number of basis shells:\t\t", length(BS1.basis))
    println("Number of basis functions: \t", lengthIndex(BS1))
    println("Length of multiindex I: \t", lengthIndex(BS1)^2,"\n")

    println("Molecule size:\n", xmax_molecule, "\n", xmin_molecule, "\n")
    println("Grid size:\n", xmax, "\n", xmin, "\n")
end

function basis()
    println("-------------------------------------------------------------------------------------")
    println("I\ti\tAtom\tx\t\ty\t\tz\t\tl\tl_xyz")
    println("-------------------------------------------------------------------------------------")
    aktuell = Vector{Float64}
    for k in 1:lengthIndex(BS1)
        i,li = getAtomL(k)
        bsf = BS1.basis[i]
        lxyz = permute(bsf.l)[li,:]
        print(k,"\t",i , "\t", bsf.atom.Z)
        if bsf.atom.xyz != aktuell
            print("\t", bsf.atom.xyz[1],"\t", bsf.atom.xyz[2],"\t", bsf.atom.xyz[3],"\t")
            aktuell = bsf.atom.xyz
        else
            print( "\t\t\t\t\t\t\t")
        end
        println(bsf.l,"\t$(lxyz[1])$(lxyz[2])$(lxyz[3])")

    end
    println("-------------------------------------------------------------------------------------")
end

"""-------------------------
Here we define our functions
-------------------------"""
#Basis function
function phi(numberBasisfunction::Int64, l_xyz::Int64, r::Vector{Float64})
    R = BS1.basis[numberBasisfunction].atom.xyz
    r = r-R

    basisfunction = BS1.basis[numberBasisfunction]
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
function a(r::Vector{Float64}, i::Int64, j::Int64)
    fakeatom = [GaussianBasis.Atom(-1, 1, r)]
    result = Coulomb_Integral(BS1, i, j, fakeatom)
    return result
end

function Coulomb_Integral(BS::BasisSet, i::Int64, j::Int64, atom::Vector{<:Atom}) #choose two Basis functions from the set 
    out = zeros(GaussianBasis.num_basis(BS.basis[i]), GaussianBasis.num_basis(BS.basis[j]))
    GaussianBasis.generate_V_pair!(out, BS.basis[i], BS.basis[j], atom)
    return out
end

#Full Discretized analytic Tensor
function A(list::Vector{Int64})
    I = list[1]
    coord = QG.quantics_to_origcoord(grid, list[2:end])

    i, li = getAtomL(splitIndex(I)[1])
    j, lj = getAtomL(splitIndex(I)[2])

    return a([coord[1], coord[2], coord[3]], i, j)[li, lj]
end

#Total Coulomb Integral (μν|λρ)
function totalIntegral(I,J)

    i,li,A1,R1,l1 = parametersPhi(splitIndex(I)[1])
    j,lj,B,R2,l2 = parametersPhi(splitIndex(I)[2])
    
    result = 0.
    for coord in Iterators.product(1:8, ntuple(x -> 1:8, R - 1)...)
        x = collect(QG.quantics_to_origcoord(grid,collect(coord)))
        y = [J,coord...]
        result += phi(A1,l1,x-R1)*phi(B,l2,x-R2)*A(y)
    end
    
    return result*prod(xmax-xmin)*ang2bohr^3/(2^(3*R))
end

"""--------------------------------------------------
Tensor cross interpolation to compress the function h
--------------------------------------------------"""

function worstcase1(L::Int64)
    worstcase = Vector{Float64}(undef,L)
    for i in 1:Int(ceil(L/2))
        worstcase[i] = 8^i
        worstcase[L+1-i] = 8^i
    end
    worstcase
end

function worstcase2(L::Int64)
    localdims = fill(8, R+1)
    localdims[1] = lengthIndex(BS1)^2

    worstcase = Vector{Float64}(undef,L)
    for i in 1:Int(L)
        worstcase[i] = min(prod(localdims[1:i]), prod(localdims[i:end]))
    end
    return worstcase
end

function fullCompress()

    localdims = fill(8, R+1)
    localdims[1] = lengthIndex(BS1)^2

    cA = TensorCrossInterpolation.CachedFunction{Float64}(A, localdims)
    tci, ranks, errors = crossinterpolate2(Float64, cA, localdims; tolerance = tol)
    r = TensorCrossInterpolation.linkdims(tci)

    return (tci,ranks,errors)

    worstcase = worstcase2(length(r))

    tick = 13
    guide = 13

    p = plot(1:length(localdims)-1, [worstcase, r], label = ["worstcase" "TCI ranks"],  yscale= :log10, legendfontsize=guide, xtickfontsize=tick, ytickfontsize=tick, xguidefontsize=guide, yguidefontsize=guide)
    title!("TT rank of Integral Function A, with R = $R")
    xlabel!("TT chain position")
    ylabel!("rank")
    savefig(p, output)
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

function atomPlot(I::Int64)
    i,li,A,R1,l1 = parametersPhi(I)
    #R1 = round.(R1, digits=2)
    xaxis1 = range(-offset, offset, length = 1000)
    xaxis2 = range(xmin[1], xmax[1], length = 1000)
    ydata1 = [phi(A,l1,[k,0.,0.]) for k in xaxis1]
    ydata2 = [phi(A,l1,[k-R1[1],0.,0.]) for k in xaxis2]
    plt1 = plot(xaxis1, ydata1, label = "Basis function $i", xlabel = "r in angstrom", ylabel = "φ(r)")
    plt2 = plot(xaxis1, ydata1, label = "Basis function $i",yscale = :log10, xlabel = "r in angstrom", ylabel = "φ(r)")
    plt3 = plot(xaxis2, ydata2, label = "Basis function $i", xlabel = "r in angstrom", ylabel = "φ(r)")
    plt4 = plot(xaxis2, ydata2, label = "Basis function $i",yscale = :log10, xlabel = "r in angstrom", ylabel = "φ(r)")
    plt = plot(plt1,plt2,plt3,plt4, layout=(2,2), legend=false, dpi = 600)
    display(plt)
    #plt = plot(plt1,plt2, layout=(1,2), legend=false, dpi = 600, title = "Basisfunction at $R1 with l_xyz = $l1",  titlefont = font(12))
end

function twoAtomPlot(I::Int64, J::Int64)
    i,li,A,R1,l1 = parametersPhi(I)
    j,lj,B,R2,l2 = parametersPhi(J)
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
function compareResults(num, ana)
    println("Numerical value:\t",num,"\tAnalytical value:\t",ana,"\t Difference:\t", num-ana)
end

function parametersPhi(n::Int64)
    i,li = getAtomL(n)
    coef1 = BS1.basis[i].coef
    exp1 = BS1.basis[i].exp
    A = hcat(coef1,exp1)
    R1 = BS1.basis[i].atom.xyz
    l1 = permute(BS1.basis[i].l)[li,:]
    return (i,li,A,R1,l1)
end

function test_A()
    I = 25 
    N = 20
    #A_tci = fullCompress()[1]

    i,li,A,R1,l1 = parametersPhi(splitIndex(I)[1])
    j,lj,B,R2,l2 = parametersPhi(splitIndex(I)[2])


    for k in 1:N

        xi = range(xmin_molecule[1], xmax_molecule[1], length = 2^R)
        yi = range(xmin_molecule[2], xmax_molecule[2], length = 2^R)
        zi = range(xmin_molecule[3], xmax_molecule[3], length = 2^R)

        x = collect(QG.quantics_to_origcoord(grid,rand(1:8,R)))
        if k ==1 
            x = collect(QG.quantics_to_origcoord(grid,fill(4,R)))
        end
        #x = [rand(xi), rand(yi), rand(zi)]

        result = 0.
        for coord in Iterators.product(1:8, ntuple(x -> 1:8, R - 1)...)
            y = collect(QG.quantics_to_origcoord(grid,collect(coord)))
            if x != y
                result += phi(A,l1,(y-R1))*phi(B,l2,(y-R2))/(norm(x-y)*ang2bohr)
            end
        end
        
        num = result*prod(xmax-xmin)*ang2bohr^3/(2^(3*R))
        ana = a(x, i,j)[li,lj]
        
        compareResults(num,ana)

    end

    #points = transpose(points)
    #plt3d = scatter(points[:,1], points[:,2], points[:,3])
    #display(plt3d)
    
end

function test_phi()
    I = 1
    N = lengthIndex(BS1)
    num = zeros(Float64, N,N)
    
    for n in 1:N
        for m in 1:N
            if m>=n 
                i,li,A,R1,l1 = parametersPhi(n)
                j,lj,B,R2,l2 = parametersPhi(m)

                #=
                result = 0.
                for ix in xs
                    for iy in ys
                        for iz in zs
                            y = [ix,iy,iz]
                            result += phi(A,l1,y-R1)*phi(B,l2,y-R2)
                        end
                    end
                end
                num[n,m] = result*norm(prod(xmax-xmin))/(2^(3*R))
                =#

                result = 0.
                for coord in Iterators.product(1:8, ntuple(x -> 1:8, R - 1)...)
                    y = collect(QG.quantics_to_origcoord(grid,collect(coord)))
                    result += phi(A,l1,(y-R1))*phi(B,l2,(y-R2))
                end
                num[n,m] = num[m,n] = result*norm(prod(xmax-xmin)*ang2bohr^3)/(2^(3*R))
            end
        end
    end

    #ana = overlap(BS1, i,j)[li,lj]
    ana = overlap(BS1, BS1)
    
    return (num-ana)
    compareResults(num,ana)
    
end

function compressPhi()

    I=1
    localdims = fill(8, R)
    N = lengthIndex(BS1)
    num = zeros(Float64, N,N)
    
    for n in 1:N
        for m in 1:N
            if m>=n 

                i,li,A,R1,l1 = parametersPhi(n)
                j,lj,B,R2,l2 = parametersPhi(m)

                p(list::Vector) = phi(A,l1,collect(QG.quantics_to_origcoord(grid, list))-R1)*phi(B,l2,collect(QG.quantics_to_origcoord(grid, list))-R2)

                cp = TensorCrossInterpolation.CachedFunction{Float64}(p, localdims)
                pivot3 = collect(QG. origcoord_to_quantics(grid,Tuple(Float64(x) for x in (R1+R2)/2)))
                tci = crossinterpolate2(Float64, cp, localdims, [pivot3]; tolerance = tol)[1]

                num[n,m] = num[m,n] = sum(tci)*norm(prod(xmax-xmin)*ang2bohr^3)/(2^(3*R))
            end
        end
    end

    ana = overlap(BS1, BS1)
    return num-ana
end

function testTotal()
    N = 20
    L = lengthIndex(BS1)^2

        
    for k in 1:N
        I,J = rand(1:L,2)

        i, li = getAtomL(splitIndex(I)[1])
        j, lj = getAtomL(splitIndex(I)[2])
        n, ln = getAtomL(splitIndex(J)[1])
        m, lm = getAtomL(splitIndex(J)[2])

        compareResults(ERI_2e4c(BS1,i,j,n,m)[li,lj,ln,lm], totalIntegral(I,J))
    end
end


