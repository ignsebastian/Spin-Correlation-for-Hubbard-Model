function initialize_state(Nup, Ndn,N)

    state = Vector{String}()
    state = ["Emp" for n in 1:N]
    p = Nup + Ndn #for iteration purpose
    nupAdditional = Nup - (Nup+Ndn-N) #for iteration purpose


    for i in N:-1:1
        if p > i
            state[i] = "UpDn"
            p -= 2
        elseif p > 0
            if i <= nupAdditional
                state[i] = "Up"
            else
                state[i] = "Dn"
            end
            p -= 1
        end
    end

    shuffle!(state)
    shuffle!(state)

    return state
end

