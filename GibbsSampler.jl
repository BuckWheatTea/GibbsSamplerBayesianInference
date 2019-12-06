using StatsBase

function read_fasta(filePath)
    file = readlines(filePath)
    seqencesCombine = string()
    for line in file
        if line == ""
            seqencesCombine*=" "
            continue
        elseif line[1] in ['A','C','G','T']
            seqencesCombine*=line
        end
    end
    return [string(seq) for seq in filter!(e->e!= "", split(seqencesCombine," "))]
end

function GammaDivide(up::Int64, down::Int64)::BigFloat
    res = BigInt(1)
    while down > up
        res*=down
        down -= 1
    end
    return 1/res
end


function GammaDivide(up::Vector{Int64}, down::Vector{Int64})::BigFloat
    # Notice: length of up and down must be equal to 4
    res = BigFloat(1)
    for i=1:4
        res*=GammaDivide(up[i], down[i])
    end
    return res
end

function getMotifs(Sequences::Vector{String}, Position::Vector{Int64}, motifLength::Int64)
    Motifs = Vector{String}()
    for i = 1:length(Sequences)
        push!(Motifs, Sequences[i][Position[i]:Position[i]+motifLength-1])
    end
    return Motifs
end


function motifScore(DnaBox::Vector{String})::Int64
    numOfDna = length(DnaBox)
    lengthOfDna = length(DnaBox[1])
    score = 0
    for i=1:lengthOfDna
        count = Dict{Char, Int64}('A'=> 0, 'C'=>0, 'G'=>0, 'T'=>0)
        for j=1:numOfDna
            count[DnaBox[j][i]] += 1
        end
        score += findmax(count)[1]
    end
    return score
end

function transposeOfStringMatrix(Sequences)
    newSequences = Vector{String}()
    numberOfSeq = length(Sequences)
    lengthOfSeq = length(Sequences[1])

    for i = 1:lengthOfSeq
        thisSeq = string()
        for j = 1:numberOfSeq
            thisSeq *= Sequences[j][i]
        end
        push!(newSequences, thisSeq)
    end
    return newSequences
end

function baseCounter(Sequence::String)
    Count = zeros(Int64, 4)
    Count[1] = count(i->i=='A',Sequence)
    Count[2] = count(i->i=='C',Sequence)
    Count[3] = count(i->i=='G',Sequence)
    Count[4] = count(i->i=='T',Sequence)
    return Count
end

function baseCounter(Sequences::Vector{String})
    Count = sum([baseCounter(sequence) for sequence in Sequences])
    return Count
end

function baseCountForEveryPosition(Sequences::Vector{String})
    Count = Vector{Vector{Int64}}()
    SeqTranspose = transposeOfStringMatrix(Sequences)
    for seq in SeqTranspose
        push!(Count, baseCounter(seq))
    end
    return Count
end

function Pi1SlideSequence(
    sequence::String,
    baseCountForAll::Vector{Int64},
    baseCountForMotifs::Vector{Int64},
    motifLength::Int64,
    alpha::Vector{Int64})

    alpha = alpha - ones(Int64, 4)
    betAllAndMotifs = baseCountForAll - baseCountForMotifs
    lengthOfthisSequence = length(sequence)
    pi1Vector = zeros(BigFloat, lengthOfthisSequence - motifLength + 1)
    for j = 1:lengthOfthisSequence - motifLength + 1
        thisMotif = sequence[j:j+motifLength-1]
        baseCountForThisMotif = baseCounter(thisMotif)
        thisPi1 = GammaDivide(betAllAndMotifs+alpha-baseCountForThisMotif, betAllAndMotifs+alpha)
        pi1Vector[j] = thisPi1
    end
    return pi1Vector
end

function Pi2SlideSequence(
    sequence::String,
    Motifs::Vector{String},
    beta::Array{Int64, 2})

    # beta = beta# - ones(Int64, 4)
    baseCountForMotifEP = baseCountForEveryPosition(Motifs)
    pi2SlideForThisSequence = ones(Int64, length(sequence) - length(Motifs[1]) + 1)
    for i=1:length(sequence) - length(Motifs[1]) + 1
        thisMotif = sequence[i:i+length(Motifs[1]) - 1]
        for j=1:length(thisMotif)
            if thisMotif[j] == 'A'
                pi2SlideForThisSequence[i]*= baseCountForMotifEP[j][1]+beta[j,:][1]
            elseif thisMotif[j] == 'C'
                pi2SlideForThisSequence[i]*= baseCountForMotifEP[j][2]+beta[j,:][2]
            elseif thisMotif[j] == 'G'
                pi2SlideForThisSequence[i]*= baseCountForMotifEP[j][3]+beta[j,:][3]
            elseif thisMotif[j] == 'T'
                pi2SlideForThisSequence[i]*= baseCountForMotifEP[j][4]+beta[j,:][4]
            end
        end
    end
    return pi2SlideForThisSequence
end

function Pi2SlideSequence(
    sequence::String,
    Motifs::Vector{String},
    beta::Array{Int64, 2},
    baseCountForAll::Vector{Int64},
    alpha::Vector{Int64}
    )

    beta = beta# - ones(Int64, 4)
    baseCountForMotifs = baseCounter(Motifs)
    baseCountForMotifsC = baseCountForAll - baseCountForMotifs
    baseCountForMotifEP = baseCountForEveryPosition(Motifs)
    pi2SlideForThisSequence = ones(Float64, length(sequence) - length(Motifs[1]) + 1)
    for i=1:length(sequence) - length(Motifs[1]) + 1
        thisMotif = sequence[i:i+length(Motifs[1]) - 1]
        for j=1:length(thisMotif)
            if thisMotif[j] == 'A'
                pi2SlideForThisSequence[i]*= (baseCountForMotifEP[j][1]+beta[j,:][1]) / (baseCountForMotifsC[1]+alpha[1])
            elseif thisMotif[j] == 'C'
                pi2SlideForThisSequence[i]*= baseCountForMotifEP[j][2]+beta[j,:][2] / (baseCountForMotifsC[2]+alpha[2])
            elseif thisMotif[j] == 'G'
                pi2SlideForThisSequence[i]*= baseCountForMotifEP[j][3]+beta[j,:][3] / (baseCountForMotifsC[3]+alpha[3])
            elseif thisMotif[j] == 'T'
                pi2SlideForThisSequence[i]*= baseCountForMotifEP[j][4]+beta[j,:][4] / (baseCountForMotifsC[4]+alpha[4])
            end
        end
    end
    return pi2SlideForThisSequence
end

function combinePi1AndPi2(Pi1Vector, Pi2Vector)
    lengthOfVec = length(Pi1Vector)
    Pi = [Pi1Vector[i]*Pi2Vector[i] for i=1:lengthOfVec]
    return Pi
end

function choosePositionByPi(Pi)
    Positions = [i for i =1:length(Pi)]
    ChoosedPos = sample(Positions, Weights(Pi))
    return ChoosedPos
end

function choosePositionByPiMax(Pi)
    # Positions = [i for i =1:length(Pi)]
    ChoosedPos = findmax(Pi)[2]
    return ChoosedPos
end

function GetMotifFromSequenceByMotifs(
    sequence::String,
    Motifs::Vector{String},
    baseCountForAll::Vector{Int64},
    alpha::Vector{Int64},
    beta::Array{Int64, 2},
    improve::Bool=false,
    ord::Bool = false
    )

    MotifLength = length(Motifs[1])
    baseCountForMotifs = baseCounter(Motifs)


    if ord
        Pi1Vec = ones(Float64, length(sequence) - MotifLength + 1)
        Pi2Vec = Pi2SlideSequence(sequence, Motifs, beta)
        Pi = combinePi1AndPi2(Pi1Vec, Pi2Vec)
        ChoosedPos = choosePositionByPi(Pi)
        motif = sequence[ChoosedPos:ChoosedPos+MotifLength-1]
        return motif
    end

    if improve
        Pi1Vec = ones(Float64, length(sequence) - MotifLength + 1)
        Pi2Vec = Pi2SlideSequence(sequence, Motifs, beta, baseCountForAll, alpha)
        Pi = combinePi1AndPi2(Pi1Vec, Pi2Vec)
        ChoosedPos = choosePositionByPi(Pi)
        motif = sequence[ChoosedPos:ChoosedPos+MotifLength-1]
        return motif
    end

    Pi1Vec = Pi1SlideSequence(
        sequence,
        baseCountForAll,
        baseCountForMotifs,
        MotifLength,
        alpha
    )
    Pi2Vec = Pi2SlideSequence(sequence, Motifs, beta)
    Pi = combinePi1AndPi2(Pi1Vec, Pi2Vec)
    ChoosedPos = choosePositionByPi(Pi)
    motif = sequence[ChoosedPos:ChoosedPos+MotifLength-1]
    return motif
end

function roundForVector(nums::Vector{Float64})
    roundN = zeros(Int64, length(nums))
    for i = 1:length(nums)
        roundN[i] = round(nums[i])
    end
    return roundN
end

function getAlphaByMotif(baseCountForAll::Vector{Int64}, baseCountForThisMotif::Vector{Int64}, numberOfSeq::Int64)
    baseCountForMotifEP = baseCountForAll - baseCountForThisMotif
    baseCountForMotifEPPer = baseCountForMotifEP / sum(baseCountForMotifEP)
    baseCountForMotifEPCD = baseCountForMotifEPPer*numberOfSeq*0.4+ones(Int64, 4)
    alpha = roundForVector(baseCountForMotifEPCD)
    return alpha
end

function getBetaByMotif(thisMotifs::Vector{String}, numberOfSeq::Int64)
    beta = zeros(Int64, length(thisMotifs[1]), 4)
    transPoseOfthisMotif = transposeOfStringMatrix(thisMotifs)
    for i = 1:length(thisMotifs[1])
        baseCountFori = baseCounter(transPoseOfthisMotif[i])
        beta[i,:] = roundForVector( baseCountFori / sum(baseCountFori) * numberOfSeq * 0.4) + ones(Int64, 4)
    end
    return beta
end

function GibbsSamplerSingleStep(
    Sequences::Vector{String},
    Motifs::Vector{String},
    baseCountForAll::Vector{Int64},
    alpha::Vector{Int64},
    beta::Array{Int64, 2},
    improve::Bool = false,
    ord::Bool = false)

    numberOfSeq = length(Sequences)
    InMotifs = copy(Motifs)
    for i = 1:numberOfSeq
        deleteat!(InMotifs, i)
        newMotif = GetMotifFromSequenceByMotifs(
            Sequences[i],
            InMotifs,
            baseCountForAll,
            alpha,beta,
            improve,
            ord
        )
        insert!(InMotifs, i, newMotif)
    end
    return InMotifs
end


function GibbsSampler(
    Sequences::Vector{String},
    MotifLength::Int64,
    iterT::Int64 = 500,
    CtrlN::Int64 = 0,
    improve::Bool = false,
    alpha::Vector{Int64} = [-1,-1,-1,-1],
    beta::Array{Int64, 2} = ones(Int64, 4, 4)*(-1),
    )

    numberOfSeq = length(Sequences)
    baseCountForAll = baseCounter(Sequences)

    # Initial alpha, beta, CtrlN if they are not given.
    if alpha[1] == -1 alpha = ones(Int64, 4)* Int64(round(numberOfSeq*0.2+1)) end
    if beta[1] == -1 beta = ones(Int64, MotifLength, 4)* Int64(round(numberOfSeq*0.2+1)) end
    if CtrlN == -1 CtrlN = floor(iterT / 2) end
    if CtrlN == 0 CtrlN = iterT+1 end
    # Initial Motif randomly
    positionChoosed = [rand(1 : length(Sequences[i])-MotifLength+1) for i=1:numberOfSeq]

    thisMotifs = [Sequences[i][positionChoosed[i]:positionChoosed[i]+MotifLength-1] for i=1:numberOfSeq]



    # recT = floor(iterT / 2)
    recT = 0
    maxMotifScore = -1
    optimalMotifs = Vector{String}()

    iterRec = Vector{Int64}()
    scoreRec = Vector{Int64}()

    for i=1:iterT
        if i < CtrlN
            thisMotifs = GibbsSamplerSingleStep(Sequences, thisMotifs, baseCountForAll, alpha, beta, improve)
        elseif i >= CtrlN
            baseCountForMotifs = baseCounter(thisMotifs)
            alpha = getAlphaByMotif(baseCountForAll, baseCountForMotifs, numberOfSeq)
            beta = getBetaByMotif(thisMotifs, numberOfSeq)
            thisMotifs = GibbsSamplerSingleStep(Sequences, thisMotifs, baseCountForAll, alpha, beta)
        end

        if i > recT
            thisMotifScore = motifScore(thisMotifs)
            if thisMotifScore > maxMotifScore
                maxMotifScore = thisMotifScore
                optimalMotifs = thisMotifs
            end
            push!(iterRec, i)
            push!(scoreRec, maxMotifScore)
        end
    end
    return optimalMotifs, maxMotifScore, iterRec, scoreRec
end

function GibbsSamplerOrdinary(
    Sequences::Vector{String},
    MotifLength::Int64,
    iterT::Int64 = 500)

    beta = ones(Int64, MotifLength, 4)
    baseCountForAll = baseCounter(Sequences)

    numberOfSeq = length(Sequences)
    positionChoosed = [rand(1 : length(Sequences[i])-MotifLength+1) for i=1:numberOfSeq]
    thisMotifs = [Sequences[i][positionChoosed[i]:positionChoosed[i]+MotifLength-1] for i=1:numberOfSeq]

    maxMotifScore = -1
    optimalMotifs = Vector{}()

    iterRec = Vector{Int64}()
    scoreRec = Vector{Int64}()


    for i=1:iterT
        # Initial Motif randomly
        thisMotifs = GibbsSamplerSingleStep(
            Sequences, thisMotifs, baseCountForAll, [2,2,2,2], beta, false, true
        )
        thisMotifScore = motifScore(thisMotifs)
        if thisMotifScore > maxMotifScore
            maxMotifScore = thisMotifScore
            optimalMotifs = thisMotifs
        end
        push!(iterRec, i)
        push!(scoreRec, maxMotifScore)
    end
    return optimalMotifs, maxMotifScore, iterRec, scoreRec
end
