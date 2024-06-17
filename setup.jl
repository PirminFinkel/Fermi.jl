using Fermi
using GaussianBasis 
using StaticArrays
using Plots
using LinearAlgebra
using QuanticsTCI
using TensorCrossInterpolation
import QuanticsGrids as QG


"""-------------------------------
Stuff we need in the next sections 
-------------------------------"""

#Create (cartesian) gaussian basis
function createBasis(c1::String, pos1::Vector{Float64}, c2::String, pos2::Vector{Float64})
    BS = BasisSet("sto-3g", """
    $c1        $(pos1[1])        $(pos1[2])        $(pos1[3])
    $c2        $(pos2[1])        $(pos2[2])        $(pos2[3])""", spherical=false)
    return BS  
end

function findPermutationOfLibrary()
    l= 3
    am1 = l
    am2 = l

    N1 = binomial(l+2,2)
    out = zeros(Int64,N1,N1)

    index1 = 1
    for l1 in am1 : -1 : 0
        for n1 in 0 : am1 - l1
            m1 = am1 - l1 - n1
            index2 = 1
            for l2 in am2 : -1 : 0
                for n2 in 0 : am2 - l2
                    m2 = am2 - l2 - n2
                    out[index1 + N1*(index2-1)] = 1000000*7+100000*l1+10000*m1+1000*n1+100*l2+10*m2+n2
                    index2 += 1
                end
            end
            index1 += 1
        end
    end
end

function permute(l::Int64)
    A = Array{Int64,2}(undef, binomial(l+2,2),3)
    i = 1
    for x in l:-1:0
        for y in l-x:-1:0
            A[i,1] = x
            A[i,2] = y
            A[i,3] = l-x-y
            i +=1
        end
    end  
    return A          
end

function getL(BS::BasisSet, l_Number::Int64)
    bsf = BS.basis[basisFunctionNumber]
    A = permute(bsf.l)
    l = A[l_Number,:]
    #var_name = string(Core.eval(Main, Symbol(basis)))
    println("Angular momentum permutation used for : $l ")
    return nothing
end

function lengthIndex(BS::BasisSet)::Int
    result = 0
    for i in 1:length(BS.basis)
        l = BS.basis[i].l
        result += factorial(l+2)/(factorial(l)*2)
    end
    return Int(result)
end

function getAtomData(i::Int64)
    a = getAtomL(i)
    basis = BS1.basis[a[1]]
    l = basis.l
    l_xyz = permute(l)[a[2],:]
    return (a[1], basis.atom.Z , l , l_xyz, basis.atom.xyz, basis.coef , basis.exp)
end 

function getAtomL(I::Int64)
    j = 0
    atom = 1

    for i in 1:length(BS1.basis)
        l = BS1.basis[i].l
        if I <= (j+factorial(l+2)/(factorial(l)*2))
            break
        else
            j += Int(factorial(l+2)/(factorial(l)*2))
            atom += 1
        end
    end
    LNumber = I-j

    return (atom,LNumber)
end

function findIndex(Rx::Float64, Ry::Float64, Rz::Float64, lx::Int64, ly::Int64, lz::Int64 )
    index = 0
    for i in 1:length(BS1.basis)
        if BS1.basis[i].atom.xyz == [Rx,Ry,Rz]
            A = permute(BS1.basis[i].l)

            for j in 1:size(A)[1]
                index += 1
                if [lx,ly,lz] == A[j,:]
                    return index
                end
            end
        else
            index += Int(factorial(BS1.basis[i].l+2)/(factorial(BS1.basis[i].l)*2))
        end
    end
    println("No Atom at this position")
    return nothing
end

function splitIndex(I::Int64)
    return (div(I-1,lengthIndex(BS1))+1, I-div(I-1,lengthIndex(BS1))*lengthIndex(BS1))
end

function convert(g::QG.DiscretizedGrid, coordinate::NTuple{N,Float64}) where {N}
    if g.includeendpoint
        all(QG.grid_min(g) .<= coordinate .<= QG.grid_max(g)) ||
            error("Bound Error: $(coordinate), min=$(grid_min(g)), max=$(QG.grid_max(g))")
    else
        all(QG.grid_min(g) .<= coordinate .< QG.grid_max(g)) ||
            error("Bound Error: $(coordinate), min=$(QG.grid_min(g)), max=$(QG.grid_max(g))")
    end
    return QG._convert_to_scalar_if_possible(
        ((coordinate .- QG.grid_min(g)) ./ QG.grid_step(g) .+ 1) .|> round .|> Int,
    )
end

function read_xyz(xyz_filename)
    out_str = ""
    natoms = 0
    open(xyz_filename) do file
        count = 1
        for line in eachline(file)
            if count==1
                natoms = parse(Int, line)
            elseif count>2 # second line is title
                out_str *= line * "\n"
            end
            count += 1
        end
    end
    return (out_str, natoms)
end

function toArray(str::String)
    daten = split(read_xyz(str)[1], "\n")
    pop!(daten)
    A = zeros(Float64, length(daten), 3)

    for i in eachindex(daten)
        v = split(daten[i], "\t" , keepempty = false)
        for j in eachindex(v)
            if j ==1 
                continue
            end
            A[i,j-1] = parse(Float64, v[j])
        end
    end
    return A
end

function splitH(str::String)
    daten = split(read_xyz(str)[1], "\n")
    natoms = read_xyz(str)[2]
    A = Array{String }(undef, natoms, 4)

    #generate Matrix
    for i in eachindex(daten)
        v = split(daten[i], " " , keepempty = false)
        for j in eachindex(v)
            A[i,j] = v[j]
        end
    end
    A = A[sortperm(A[:, 1]), :]
 
    #Split the Matrix
    h = 0
    for i in eachindex(A[:,1])
        if A[i,1] == "H" 
            h = i
            break
        end
    end

    o = 0
    for i in eachindex(A[:,1])
        if A[i,1] == "O" 
            o = i
            break
        end
    end
    
    #Change Type to Float64
    A = A[:,2:end]
    
    B = Array{Float64}(undef,size(A))
    for i in eachindex(A)
        B[i] = parse(Float64, A[i])
    end

    C = B[1:h-1,:] 
    H = B[h:o-1,:] 
    O = B[o:end,:] 

    return (C,H,O)
end

function getMaxima(A::Array)
    vmax = [maximum(A[:, i]) for i in 1:3]
    vmin = [minimum(A[:, i]) for i in 1:3]
    return (vmax,vmin)
end

function generateOutputName(input_name::String, i::Int64)
    return s = input_name[1:end-4]*"_$i"*".png"
end

function determineGridBounds(str::String)
    xmin = getMaxima(toArray(input))[2]
    xmax = getMaxima(toArray(input))[1]
    a = abs.(xmax-xmin)

end