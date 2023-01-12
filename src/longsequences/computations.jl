
###
### Computations
###
###
### Computation methods for LongSequence
###
### This file is a part of BioJulia.
### License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

"""
    molecular_weight(seq::LongRNA, double_stranded::Bool, circular::Bool, Monoisotopic::Bool)

Calulates the molecular weight of a LongSequence
"""
###
### Molecular weight for RNA
###

function molecular_weight(
    seq::LongRNA;
    double_stranded::Bool=false,
    circular::Bool=false,
    monoisotopic::Bool=false)

    weight = 0

    if monoisotopic
        weight_table = monoisotopic_unambiguous_rna_weights
        water = 18.010565
    else
        weight_table = unambiguous_rna_weights
        water = 18.0153
    end

    try
        weight = sum(weight_table[x] for x in seq) - (length(seq) - 1) * water
        if circular
            weight -= water
        end
    catch e
        println("There is a not valid unambiguos letter for $seq_type")
    end


    if double_stranded
        seq = complement(seq)
        weight += sum(weight_table[x] for x in seq) - (length(seq) - 1) * water
        if circular
            weight -= water
        end
    end

    return weight

end

###
### Molecular weight for DNA
###
function molecular_weight(
    seq::LongDNA;
    double_stranded::Bool=false,
    circular::Bool=false,
    monoisotopic::Bool=false)

    weight = 0

    if monoisotopic
        weight_table = monoisotopic_unambiguous_dna_weights
        water = 18.010565
    else
        weight_table = unambiguous_dna_weights
        water = 18.0153
    end

    try
        weight = sum(weight_table[x] for x in seq) - (length(seq) - 1) * water
        if circular
            weight -= water
        end
    catch e
        println("There is a not valid unambiguos letter for $seq_type")
    end


    if double_stranded
        seq = complement(seq)
        weight += sum(weight_table[x] for x in seq) - (length(seq) - 1) * water
        if circular
            weight -= water
        end
    end

    return weight

end

###
### Molecular weight for AA
###
function molecular_weight(
    seq::LongAA;
    double_stranded::Bool=false,
    circular::Bool=false,
    monoisotopic::Bool=false)

    weight = 0

    if monoisotopic
        weight_table = monoisotopic_protein_weights
        water = 18.010565
    else
        weight_table = protein_weights
        water = 18.0153
    end

    try
        weight = sum(weight_table[x] for x in seq) - (length(seq) - 1) * water
        if circular
            weight -= water
        end
    catch e
        println("There is a not valid unambiguos letter for $seq_type")
    end


    if double_stranded
        throw(ErrorException("protein sequences cannot be double-stranded"))
    end

    return weight

end


# Average masses of monophosphate deoxy nucleotides
unambiguous_dna_weights = Dict{DNA, Number}(
    DNA_A=> 331.2218,
    DNA_C=> 307.1971,
    DNA_G=> 347.2212,
    DNA_T=> 322.2085
)

# Monoisotopic masses of monophospate deoxy nucleotides
monoisotopic_unambiguous_dna_weights = Dict{DNA, Number}(
    DNA_A=> 331.06817,
    DNA_C=> 307.056936,
    DNA_G=> 347.063084,
    DNA_T=> 322.056602,
)

# Average masses of monophospate nucleotides
unambiguous_rna_weights = Dict{RNA, Number}(
    RNA_A=> 347.2212,
    RNA_C=> 323.1965,
    RNA_G=> 363.2206,
    RNA_U=> 324.1813
)

# Monoisotopic masses of monophospate nucleotides
monoisotopic_unambiguous_rna_weights = Dict{RNA, Number}(
    RNA_A=> 347.063084,
    RNA_C=> 323.051851,
    RNA_G=> 363.057999,
    RNA_U=> 324.035867,
)

# Average masses of amino acids
protein_weights = Dict{AminoAcid, Number}(
    AA_A=> 89.0932,
    AA_C=> 121.1582,
    AA_D=> 133.1027,
    AA_E=> 147.1293,
    AA_F=> 165.1891,
    AA_G=> 75.0666,
    AA_H=> 155.1546,
    AA_I=> 131.1729,
    AA_K=> 146.1876,
    AA_L=> 131.1729,
    AA_M=> 149.2113,
    AA_N=> 132.1179,
    AA_O=> 255.3134,
    AA_P=> 115.1305,
    AA_Q=> 146.1445,
    AA_R=> 174.201,
    AA_S=> 105.0926,
    AA_T=> 119.1192,
    AA_U=> 168.0532,
    AA_V=> 117.1463,
    AA_W=> 204.2252,
    AA_Y=> 181.1885,
)

# Monoisotopic masses of amino acids
monoisotopic_protein_weights = Dict{AminoAcid, Number}(
    AA_A=> 89.047678,
    AA_C=> 121.019749,
    AA_D=> 133.037508,
    AA_E=> 147.053158,
    AA_F=> 165.078979,
    AA_G=> 75.032028,
    AA_H=> 155.069477,
    AA_I=> 131.094629,
    AA_K=> 146.105528,
    AA_L=> 131.094629,
    AA_M=> 149.051049,
    AA_N=> 132.053492,
    AA_O=> 255.158292,
    AA_P=> 115.063329,
    AA_Q=> 146.069142,
    AA_R=> 174.111676,
    AA_S=> 105.042593,
    AA_T=> 119.058243,
    AA_U=> 168.964203,
    AA_V=> 117.078979,
    AA_W=> 204.089878,
    AA_Y=> 181.073893,
)