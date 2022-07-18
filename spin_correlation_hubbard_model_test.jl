#Spin Correlation for Hubbard Model
#Julia version of psharma spin correlation for hubbard model
#triangular lattice

#------------------------------------Importing packages--------------------------------
using ITensors
using Random
using LinearAlgebra
using BenchmarkTools

include("bondPairs.jl")
include("initializing_state.jl")

Random.seed!(1123)

#------------------------------------------Set multithreading--------------------------
ITensors.Strided.set_num_threads(1)
BLAS.set_num_threads(1)
ITensors.enable_threaded_blocksparse()

    @show Threads.nthreads()
    @show Sys.CPU_THREADS
    @show ITensors.blas_get_num_threads()
    @show ITensors.Strided.get_num_threads()
    @show ITensors.using_threaded_blocksparse()


@time let
    #-----------------------------------------Parameter------------------------------------
    Nx = 6
    Ny = 4
    N = Nx*Ny

    Nup = 12
    Ndn = 12

    t = 1
    U = 20
    h = 0.2

    bond_dim = 10 #Bond dimension for initial state

    #-----------------------Preparing Hamiltonian on the triangular lattice-------------------#

    bonds = bondPairs(Nx,Ny) #Defining nearest neighbor pairs

    sites = siteinds("Electron", N; conserve_qns=true)

    println("Triangular XC6 lattice NN pair")

    for bnd in bonds
        println("Bond from side ", bnd.first, " -> ", bnd.second)
    end

    ampo = OpSum()
    for bnd in bonds
        ampo += -t,"Cdagup",bnd.first,"Cup",bnd.second
        ampo += -t,"Cdagup",bnd.second,"Cup",bnd.first
        ampo += -t,"Cdagdn",bnd.first,"Cdn",bnd.second
        ampo += -t,"Cdagdn",bnd.second,"Cdn",bnd.first
    end

    for i in 1:N
        ampo += U,"Nupdn",i
    end

    H = MPO(ampo,sites)
    H = splitblocks(linkinds,H)

    #= add S2 here
    for i in range(N)
        ampo += -h,"S2",i
    end

    HS2 = MPO(ampo,sites)
    =#
    #------------------------Initializing wavefunction MPS------------------------------------#

    state = initialize_state(Nup,Ndn,N)
    psi0 = randomMPS(sites,state;linkdims=bond_dim)

    #------------------------------------------------------------------------------------------
    #Setting sweep parameters
    #------------------------------------------------------------------------------------------

    #=
    println("Warmup sweeps\n")
    sweeps = Sweeps(11)
    setmaxdim!(sweeps, 8,16,24,50,100,250,500,1000,1000,2000,2000)
    setcutoff!(sweeps, 1E-8)
    noise!(sweeps, 1E-7,1E-8,0.0)
    println(sweeps)

    #Warming up initial state

    #energy, psi1 = dmrg(HS2,psi0,fsweeps,outputlevel=0; eigsolve_maxiter = 4)
    #println("Done warming up\n")

    #Final sweeps
    #
    fsweeps = Sweeps(20)
    setmaxdim!(fsweeps, 2000,3000,3000,3000,5000,5000,5000,7000,7000,7000,7000,10000,10000,10000,10000,10000,140000)
    setcutoff!(fsweeps, 1E-8)
    noise!(fsweeps, 1E-7, 1E-8,0.0)
    =#

    sweeps = Sweeps(35)
    setmaxdim!(sweeps, 8,16,24,50,100,500,1000,1000,2000,2000,3000,3000,3000,5000,5000,5000,7000,7000,7000,10000,10000,10000,10000,14000,14000,14000,14000,18000,18000,18000,18000,18000,18000,18000,18000)
    setcutoff!(sweeps, 1E-8)
    noise!(sweeps, 1E-7,1E-8,0.0)
    println(sweeps)

    #------------------------------------------------------------------------------------------
    #DMRG calculation for ground state
    #------------------------------------------------------------------------------------------

    energy, psi = @time dmrg(H,psi0,sweeps,outputlevel=1; eigsolve_maxiter=4)

    println("Initial energy = ", inner(psi0',H,psi0))
    println("\nGround State Energy from DMRG= ", energy)
    println("\nEnergy using overlap = ", inner(psi',H,psi))
    println("\nEnergy per site = ",energy/N)

    E2 = inner(H,psi,H,psi)

    println("\n<ψ|H^2|ψ> = ",E2)
    println("\nvariance = ", E2-energy*energy)
    println("\noverlap with initial state = ", inner(psi0,psi))
    println("\nTotal QN of Ground State = ", totalqn(psi))

    #------------------------------------------------------------------------------------------
    #Saving MPS to a file
    #------------------------------------------------------------------------------------------


    #------------------------------------------------------------------------------------------
    #Calculating <Nup_i> <Ndn_i> <Nupdn_i>
    #------------------------------------------------------------------------------------------

    Cup = expect(psi,"Nup")
    Cdn = expect(psi,"Ndn")
    Cupdn = expect(psi,"Nupdn")

    col_name = "i \t <Nup_i> \t <Ndn_i> \t <Nupdn_i> \n"
    println(col_name)

    output_density_corr = "density_correlator.txt"
    open(output_density_corr,"w") do file
    write(file, col_name)
        for i in 1:length(Cup)
            data_text = "$i \t $(Cup[i]) \t $(Cdn[i]) \t $(Cupdn[i])\n"
            println(" ",i," ",Cup[i]," ",Cdn[i]," ",Cupdn[i])
            write(file, data_text)
        end
    end

    #------------------------------------------------------------------------------------------
    #Calculating magnetization
    #------------------------------------------------------------------------------------------
    magz = expect(psi, "Sz")
    @show magz
    output_magz = "magnetization.txt"
    open(output_magz,"w") do file
        for i in 1:size(magz)[1]
            text = "$i \t $(magz[i])\n"
            write(file, text)
        end
    end

    #------------------------------------------------------------------------------------------
    #Calculating spin correlator
    #------------------------------------------------------------------------------------------

    C = correlation_matrix(psi,"Sz","Sz")
    Cpm = correlation_matrix(psi,"S+","S-")
    Cmp = correlation_matrix(psi,"S-","S+")

    col_name = "i \t j \t <SzSz> \t <S+S-> \t <S-S+> \n"
    println(col_name)
    output_spin_corr = "spin_correlator.txt"
    open(output_spin_corr,"w") do file
        write(file,col_name)
        for i in 1:size(C)[1]
            for j in i:size(C)[1]
                data_text = "$i \t $j \t $(C[i,j]) \t $(Cpm[i,j]) \t $(Cmp[i,j])\n"
                println(" ",i," \t",j," \t",C[i,j]," \t",Cpm[i,j]," \t",Cmp[i,j])
                write(file, data_text)
            end
        end
    end

end
