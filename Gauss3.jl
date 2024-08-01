include("setup.jl")
include("nicePlots.jl")

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

    result = 0.
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
function a(BS::BasisSet, r::AbstractVector, i::Int64, j::Int64)
    fakeatom = [GaussianBasis.Atom(-1, 1, r)]
    result = Coulomb_Integral(BS, i, j, fakeatom)
    return result
end

function Coulomb_Integral(BS::BasisSet, i::Int64, j::Int64, atom::Vector{<:GaussianBasis.Atom}) #choose two Basis functions from the set 
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

function generalKernel(c::Config, list::Vector{Int64})::Float64
    x = collect(QG.quantics_to_origcoord(c.grid,collect(list[1:c.grid.R])))
    y = collect(QG.quantics_to_origcoord(c.grid,collect(list[c.grid.R+1:2*c.grid.R])))
    return u(x,y)
end

function u( x::Vector{Float64}, y::Vector{Float64})
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

function fullCompress(c::Config, tol::Float64)
    R = c.grid.R
    localdims = fill(8, R+1)
    L = c.BS.nbas
    localdims[1] = L^2

    A_comp(list::Vector) = A(c, list)

    #R1 = Tuple(c.BS.basis[getAtomL(c.BS, L)[1]].atom.xyz)
    #pivot = [collect(QG.origcoord_to_quantics(c.grid, R1))]

    cA = TCI.CachedFunction{Float64}(A_comp, localdims)
    tci, ranks, errors = TCI.crossinterpolate2(Float64, cA, localdims; tolerance = tol) #initialpivot [[1,1,...]] without any keyword
    return TCI.linkdims(tci)
end

function makeTT(c::Config, tol::Float64)
    localdims = fill(8, c.grid.R+1)
    localdims[1] = c.BS.nbas^2

    A_comp(list::Vector{Int64}) = A(c, list) 
    P_comp(list::Vector{Int64}) = P(c, list) 
    
    #Tensor A
    cA = TCI.CachedFunction{Float64}(A_comp, localdims)
    tci_A = crossinterpolate2(Float64, cA, localdims; tolerance = tol)[1]

    #Tensor P 
    cP = TCI.CachedFunction{Float64}(P_comp, localdims)
    tci_P = crossinterpolate2(Float64, cP, localdims; tolerance = tol)[1]

    return (tci_A, tci_P)
end

function makeGeneralTT(c::Config, P_general::Function, K_general::Function, basislength::Int64)
    R = c.grid.R
    localdims = fill(8, R+1)
    localdims[1] = basislength^2
    kerneldims = fill(8,2*R)

    #Tensor P
    P_comp(list::Vector{Int64}) = P_general(c, list)
    cP = TCI.CachedFunction{Float64}(P_comp, localdims)
    tci_P, ranks_P, errors_P = crossinterpolate2(Float64, cP, localdims; tolerance = tol)

    #Tensor K 
    K_comp(list::Vector{Int64}) = K_general(c, list)
    cK = TCI.CachedFunction{Float64}(K_comp, kerneldims)
    tci_K, ranks_K, errors_K = crossinterpolate2(Float64, cK, kerneldims; tolerance = tol)

    return (tci_P, tci_K)
end

function evaluateTT(grid::QG.Grid, tci_A, tci_P)
    R = grid.R

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

    return Matrix(result, inds(result)[1], inds(result)[2])*prod(collect(QG.grid_max(grid))-collect(QG.grid_min(grid)))*ang2bohr^3/(2^(3*R))
    #return Matrix(M, inds(M)[1], inds(M)[2])*norm(prod(xmax-xmin)*ang2bohr^3)/(2^(3*R))
end

function evaluateGeneralTT(c::Config, tci_A, tci_B, tci_C)
    R = c.grid.R
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
    return Matrix(result, inds(result)[1], inds(result)[2])*prod(collect(QG.grid_max(grid))-collect(QG.grid_min(grid)))*ang2bohr^3/(2^(3*R))
end

function compress_P(c::Config, tol::Float64)

    R = c.grid.R

    localdims = fill(8,R+1)
    localdims[1] = c.BS.nbas^2
    fun(list::Vector{Int64}) = P(c, list)

    firstpivot = ones(Int,R+1)
    firstpivot[2] = 8
    TCI.optfirstpivot(fun, localdims, firstpivot )

    cf = TCI.CachedFunction{Float64}(fun, localdims)
    tci, ranks, errors = TCI.crossinterpolate2(Float64, cf, localdims, [firstpivot]; tolerance = tol) #initialpivot [[1,1,...]] without any keyword
    #print("Ranks: ", TCI.linkdims(tci))
    return tci
end

function compress_P_fixed_orbital(c::Config, I::Int64, tol::Float)
    R = c.grid.R

    localdims = fill(8,R)
    fun(list::Vector{Int64}) = P(c,[ I, list...])

    firstpivot = ones(Int,R)
    firstpivot[1] = 8
    TCI.optfirstpivot(fun, localdims, firstpivot )

    cf = TCI.CachedFunction{Float64}(fun, localdims)
    tci, ranks, errors = TCI.crossinterpolate2(Float64, cf, localdims, [firstpivot]; tolerance = tol) #initialpivot [[1,1,...]] without any keyword
    #print("Ranks: ", TCI.linkdims(tci))
    return tci
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
        println(v[k])

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

                cp = TCI.CachedFunction{Float64}(p, localdims)
                pivot3 = collect(QG. origcoord_to_quantics(c.grid,Tuple(Float64(x) for x in (R1+R2)/2)))
                tci = crossinterpolate2(Float64, cp, localdims, [pivot3]; tolerance = tol)[1]

                num[n,m] = num[m,n] = sum(tci)*prod(collect(QG.grid_max(c.grid))-collect(QG.grid_min(c.grid)))*ang2bohr^3/(2^(3*R))
            end
        end
    end

    ana = overlap(c.BS, c.BS)
    return num-ana
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


function largeMolecules(R::Int64)
    mol = ("hydrogen","water","ethane", "propanol", "glucose")
    #mol = ("hydrogen","ethane", "propanol", "glucose")
    #mol = ("hydrogen",)
    v = Vector{Config}(undef,2*length(mol))

    println("Create basis sets and grids (PF)")
    for n in eachindex(mol)
        for m in (1,2)
            v[2*(n-1)+m] =  Config(mol[n]*".txt",m, R)
            #v[4*(n-1)+m+2] = Config(mol[n]*".txt",m,7)
        end
    end

    daten = Dict() 
    println("Run has started (PF)")
    for n in eachindex(v)
        c = v[n]
        time = @elapsed ranks = fullCompress(c,1e-5)
        daten[c.name*"$(c.grid.R)"*"$(2-n%2)"] = (time,ranks,c.BS.nbas)
        println(c.name*"$(c.grid.R)"*"$(2-n%2)  with the runtime: $(time) s and the ranks: ", ranks, "Number of Orbitals: ", c.BS.nbas); flush(stdout)
        #run(`echo "$(c.name)$(c.grid.R)$(2-n%2) with the runtime: $(time) s and the ranks: $(ranks)"`)
    end

    open("daten.json","w") do file
        write(file, JSON.json(daten))
    end

    println("Run successfull (PF)")
end

function compress_multiple_P(name::String, lib::Int64, tol::Int64)
    
    daten = Dict()
    for R in (4,6,8)
        c = Config(name, lib, R)
        time = @elapsed ranks = TCI.linkdims(compress_P(c, 10.0^tol)) 
        daten["$R"] = (time, ranks)
        println("$(c.name), lib = $lib, nbas = $(c.BS.nbas), R = $(c.grid.R), tol = $tol, time = $time, ranks =\t$ranks"); flush(stdout) 
    end

    open("./daten/P_$(name)$(tol).json","w") do file
        write(file, JSON.json(daten))
    end

    return daten
end
 
function myplot()
    c = Config("water.txt",2,5)

    #plt = plot_function_axis(c.BS, 1)
    plt = plot_P(3,16,c.BS) 
    #plt = plot_multiple_in_x(c.BS, (2,3,4,7,10), "O")
    #plt = plot_multiple_in_x(c.BS, (16,17,18), "H")
    display(plt)
end

function main()
    println("Start programm") 

    compress_multiple_P("propanol.txt", 2, -5)

    println("run successfull")
end
 
#main()
myplot()
#largeMolecules(5)
#test_json()
