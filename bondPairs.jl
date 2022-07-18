#Creating bond pairs for triangular lattice
#Julia version of psharma code

function bondPairs(Nx, Ny)
    N = Nx*Ny
    bonds = Vector{Pair}()

    for i in 1:N-1
        if mod(i,Ny) == 1 && i<N-Ny
            p1 = Pair(i,i+1)
            push!(bonds,p1)
            p1 = Pair(i,i+Ny)
            push!(bonds,p1)
            p1 = Pair(i,i+Ny+1)
            push!(bonds,p1)
            p1 = Pair(i,i+Ny-1)
            push!(bonds,p1)
            p1 = Pair(i,i+2*Ny-1)
            push!(bonds,p1)

        elseif mod(i,Ny) == 0 && i <N
            p1 = Pair(i,i+Ny)
            push!(bonds,p1)

        elseif i>N-Ny
            p1 = Pair(i,i+1)
            push!(bonds,p1)
            if i==N-Ny+1
                p1 = Pair(i,i+Ny-1)
                push!(bonds,p1)
            end

        elseif mod(i,2) == 1
            p1 = Pair(i,i+1)
            push!(bonds,p1)
            p1 = Pair(i,i+Ny-1)
            push!(bonds,p1)
            p1 = Pair(i,i+Ny)
            push!(bonds,p1)
            p1 = Pair(i,i+Ny+1)
            push!(bonds,p1)

        elseif mod(i,2) == 0
            p1 = Pair(i,i+1)
            push!(bonds,p1)
            p1 = Pair(i,i+Ny)
            push!(bonds,p1)

        end
    end

    bonds = union(bonds) #to remove duplicates

    return bonds

end
