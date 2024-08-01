# using Fermi
using GaussianBasis 
using StaticArrays
using LinearAlgebra
using QuanticsTCI
using JSON

import QuanticsGrids as QG
import TensorCrossInterpolation as TCI

#=
using Profile
using StatProfilerHTML
using TCIITensorConversion
using ITensors
=#

"""-----------------------------------------------------------------------------------------------
Helper functions
-----------------------------------------------------------------------------------------------"""

function parametersPhi(BS::BasisSet, n::Int64)::Tuple{Int64, Int64, Matrix, Vector, Vector{Int64}}
    i,li = getAtomL(BS, n)
    coef1 = BS.basis[i].coef
    exp1 = BS.basis[i].exp
    A = hcat(coef1,exp1)
    R1 = BS.basis[i].atom.xyz
    l1 = permute(BS.basis[i].l)[li,:]
    return (i,li,A,R1,l1)
end

function compareResults(num::Float64, ana::Float64)
    println("Numerical value:\t",num,"\tAnalytical value:\t",ana,"\t Difference:\t", num-ana)
end

function compareResults2(avg::Float64, max::Float64)
   println("Average diviation:\t", avg, "\nMaximum diviation:\t", max) 
end

function generateOutputName(input_name::String, R::Int64, tol)
    return input_name*"_$(R)_$(tol)"*".png"
end

"""-----------------------------------------------------------------------------------------------
Permutation
-----------------------------------------------------------------------------------------------"""

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

"""-----------------------------------------------------------------------------------------------
Index 
-----------------------------------------------------------------------------------------------"""

function getL(BS::BasisSet, l_Number::Int64)
    bsf = BS.basis[basisFunctionNumber]
    A = permute(bsf.l)
    l = A[l_Number,:]
    #var_name = string(Core.eval(Main, Symbol(basis)))
    println("Angular momentum permutation used for : $l ")
    return nothing
end

function lengthIndex(BS::BasisSet)::Int
    return BS.nbas
end

function getAtomL(BS::BasisSet, I::Int64)
    j = 0
    atom = 1

    for i in 1:length(BS.basis)
        l = BS.basis[i].l
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

function splitIndex(BS::BasisSet, I::Int64)
    return (div(I-1,lengthIndex(BS))+1, I-div(I-1,lengthIndex(BS))*lengthIndex(BS))
end

"""-----------------------------------------------------------------------------------------------
Read input files and convert strings to arrays
-----------------------------------------------------------------------------------------------"""

function read_json(filename::String)
    json_string = read(filename, String)
    return JSON.parse(json_string)
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


"""-----------------------------------------------------------------------------------------------
Grid and bounds
-----------------------------------------------------------------------------------------------"""

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

function getMaxima(A::Array)
    vmax = [maximum(A[:, i]) for i in 1:3]
    vmin = [minimum(A[:, i]) for i in 1:3]
    return (vmax,vmin)
end

function getOffset(BS::BasisSet, min::Vector{Float64}, max::Vector{Float64}, limit::Float64)

    incr = 0.01

    Off = Matrix{Float64}(undef, 2,3)
    start = Vector{Float64}

    for n in (-3,-2,-1,1,2,3)
        dir = [0.,0.,0.]   
        dir[abs(n)] = sign(n)
        if n<0
            start = min
        else
            start = max
        end


        offset = 0.0
        enough = false
        while !enough
            offset += incr
            enough = true 
            for m in 1:lengthIndex(BS)
                i,li,A,R1,l1 = parametersPhi(BS, m)

                if abs(phi(A,l1,(start[abs(n)]-R1[abs(n)]+sign(n)*offset)*dir)) > limit
                    enough = false
                    break
                end
            end

            if enough
                break
            end
        end
        
        Off[Int(1.5+0.5*sign(n)), abs(n)] = offset
    end 

    Off = round.(Off, digits =2)
    offsetmin = Off[1,:]
    offsetmax = Off[2,:]
    return (offsetmin, offsetmax)
end

function get_perfect_P_grid(I::Int64, J::Int64, BS::BasisSet)
     
    i,li,A,R1,l1 = parametersPhi(BS, J)
    j,lj,B,R2,l2 = parametersPhi(BS, I)
    f1(r) = phi(A, l1, r)
    f2(r) = phi(B, l2, r)

    incr = 0.1

    Off = zeros(Float64, 2,3)
    start = Vector{Float64}

    for n in (0,1,2)
        for m in (1,2)
            sig = 2*m-3

            offset = 0.0
            enough = false
            while !enough
                Off[n%3+1,m] += incr*sig
                enough = true 
                
                for k1 in range(Off[(n+1)%3+1,1], Off[(n+1)%3+1,2], length = 5)
                    for k2 in range(-Off[(n+2)%3+1,1], Off[(n+2)%3+1,2], length = 5)

                        r = [0.,0.,0.]
                        r[n%3+1] = Off[n%3+1,m]
                        r[(n+1)%3+1] = k1 
                        r[(n+2)%3+1] = k2 

                        if abs(P(r)) > limit
                            enough = false
                            break
                        end
                        if !enough
                            break
                        end
                    end
                    if !enough
                        break 
                    end
                end

                if enough
                    break
                end
            end
            
        end 
    end

    return Off
    offsetmin = Off[:,1]
    offsetmax = Off[:,2]
    return (offsetmin, offsetmax)
    
end

function grid(BS::BasisSet, input::String, R::Int64)
    tol = 10e-7
    xmin_molecule = getMaxima(toArray(input))[1]
    xmax_molecule = getMaxima(toArray(input))[2]
    xmin = xmin_molecule - getOffset(BS, xmin_molecule, xmax_molecule, tol)[1]
    xmax = xmax_molecule + getOffset(BS, xmin_molecule, xmax_molecule, tol)[2]

    #shift = (xmax-xmin)/((2^R)-1)/(10^6)
    return QG.DiscretizedGrid{3}(R, (xmin[1], xmin[2], xmin[3]), (xmax[1], xmax[2], xmax[3]); includeendpoint=true)
end


"""-----------------------------------------------------------------------------------------------
Class Config and constructors for grid and Basis 
-----------------------------------------------------------------------------------------------"""

struct Config
    name::String
    BS::BasisSet
    grid::QG.Grid
end

function BasisSetConstructor(input::String, lib::Int64)
    whichBasis = ("sto-3g", "def2-svp") 
    return BasisSet(whichBasis[lib], read_xyz(input)[1], spherical=false, lib =:acsint)
end

function Config(input::String, lib::Int64, R::Int64)
    name = input[1:end-4]
    BS = BasisSetConstructor(input,lib)
    g = grid(BS, input, R)
    Config(name, BS, g)
end
