function simulateSequence(length::Int64, distribution::Vector{Float64} = [0.25,0.25,0.25,0.25])
    sequence = string()
    for i=1:length
        sampleBase = sample(['A','C','G','T'], Weights(distribution))
        sequence*=sampleBase
    end
    return sequence
end

function simulateMultiSequences(lengths::Vector{Int64}, distribution::Vector{Float64} = [0.25,0.25,0.25,0.25])
    Sequences = Vector{String}()
    for length in lengths
        push!(Sequences, simulateSequence(length, distribution))
    end
    return Sequences
end

function InsertMotifsToSequences(Sequences::Vector{String}, Motif::String)
    newSequences = Vector{String}()
    Motifs = Vector{String}()
    for i = 1:length(Sequences)
        thisMotif = Motif
        lengthOfthisSequence = length(Sequences[i])
        pos = rand(0:lengthOfthisSequence+1)
        if 0 < pos <  lengthOfthisSequence+1
            newSequence = Sequences[i][1:pos]*thisMotif*Sequences[i][pos+1:lengthOfthisSequence]
        elseif pos == 0
            newSequence = thisMotif*Sequences[i]
        elseif pos == lengthOfthisSequence+1
            newSequence = Sequences[i]*thisMotif
        end
        push!(newSequences, newSequence)
        push!(Motifs, thisMotif)
    end
    return newSequences, Motifs
end

function InsertMotifsToSequences(Sequences::Vector{String}, Motifs::Vector{String})
    newSequences = Vector{String}()
    for i = 1:length(Sequences)
        thisMotif = Motifs[i]
        lengthOfthisSequence = length(Sequences[i])
        pos = rand(0:lengthOfthisSequence+1)
        if 0 < pos <  lengthOfthisSequence+1
            newSequence = Sequences[i][1:pos]*thisMotif*Sequences[i][pos+1:lengthOfthisSequence]
        elseif pos == 0
            newSequence = thisMotif*Sequences[i]
        elseif pos == lengthOfthisSequence+1
            newSequence = Sequences[i]*thisMotif
        end
        push!(newSequences, newSequence)
    end
    return newSequences, Motifs
end


function OneOnOneCompare(Motif1::String, Motif2::String)
    MotifLength = length(Motif1)
    cnt = 0
    for i = 1:MotifLength
        cnt += Int64(Motif1[i] == Motif2[i])
    end
    cnt /= MotifLength
    return cnt
end

function QuietSameCompare(Motif1::String, Motif2::String)
    return Int64(Motif1 == Motif2)
end


function globalAlignmentCompare(Motif1::String, Motif2::String)
    matchScore = 1
    mismatchScore = -1
    indelScore = 0

    lengthOfMotif = length(Motif1)
    nMatrix = zeros(Int64, lengthOfMotif+1, lengthOfMotif+1)
    for i = 2:lengthOfMotif+1
        nMatrix[i, 1] = nMatrix[i-1,1] + indelScore
    end
    for j = 2:lengthOfMotif+1
        nMatrix[1, j] = nMatrix[1,j-1] + indelScore
    end
    nMatrix[1,1] = 0
    for i=1:lengthOfMotif
        for j = 1:lengthOfMotif
            if Motif1[i] == Motif2[j]
                match = matchScore
            elseif Motif1[i] != Motif2[j]
                match = mismatchScore
            end
            nMatrix[i+1,j+1] = max(nMatrix[i+1,j] + indelScore, nMatrix[i,j+1]+ indelScore, nMatrix[i,j] + match)
        end
    end
    return nMatrix[lengthOfMotif+1,lengthOfMotif+1] / lengthOfMotif
end



function CompareMotifs(Compare, Motifs1::Vector{String}, Motifs2::Vector{String})
    numOfMotif = length(Motifs1)
    score = Float64(0)
    for i = 1:numOfMotif
        score += Compare(Motifs1[i], Motifs2[i])
    end
    score = score / numOfMotif
    return score
end
