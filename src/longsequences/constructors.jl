###
### Constructors
###
###
### Constructor methods for LongSequences.
###
### This file is a part of BioJulia.
### License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

function LongSequence{A}(len::Integer) where {A<:Alphabet}
    if len < 0
        throw(ArgumentError("len must be non-negative"))
    end
    return LongSequence{A}(Vector{UInt64}(undef, seq_data_len(A, len)), convert(Int, len))
end

LongSequence(::Type{DNA}) = LongDNASeq()
LongSequence(::Type{RNA}) = LongRNASeq()
LongSequence(::Type{AminoAcid}) = LongAminoAcidSeq()
LongSequence(::Type{Char}) = LongCharSeq()

function LongSequence()
    return LongSequence{VoidAlphabet}(Vector{UInt64}(), 0)
end

function LongSequence{A}(seq::LongSequence{A}, part::UnitRange) where A
    return seq[part]
end

function LongSequence{A}(s::Union{String, SubString{String}}) where {A<:Alphabet}
    return LongSequence{A}(s, codetype(A()))
end

# Generic method for String/Substring.
function LongSequence{A}(s::Union{String, SubString{String}}, ::AlphabetCode) where {A<:Alphabet}
    len = length(s)
    seq = LongSequence{A}(len)
    return copyto!(seq, 1, s, 1, len)
end

function LongSequence{A}(s::Union{String, SubString{String}}, ::AsciiAlphabet) where {A<:Alphabet}
    v = GC.@preserve s unsafe_wrap(Vector{UInt8}, pointer(s), ncodeunits(s))
    seq = LongSequence{A}(length(v))
    return encode_chunks!(seq, 1, v, 1, length(v))
end

function LongSequence{A}(
        src::Union{AbstractString,AbstractVector},
        startpos::Integer=1,
        stoppos::Integer=length(src)) where {A<:Alphabet}
    len = stoppos - startpos + 1
    seq = LongSequence{A}(len)
    return copyto!(seq, 1, src, startpos, len)
end

# create a subsequence
function LongSequence(other::LongSequence, part::UnitRange{<:Integer})
    checkbounds(other, part)
    subseq = typeof(other)(length(part))
    copyto!(subseq, 1, other, first(part), length(part))
    return subseq
end

# Create a 4 bit DNA/RNA sequences from a 2 bit DNA/RNA sequences, and vice-versa.
for (alpha, alphb) in [(DNAAlphabet{4}, DNAAlphabet{2}), # DNA to DNA
                       (DNAAlphabet{2}, DNAAlphabet{4}),
                       (RNAAlphabet{4}, RNAAlphabet{2}), # RNA to RNA
                       (RNAAlphabet{2}, RNAAlphabet{4}),
                       (DNAAlphabet{2}, RNAAlphabet{4}), # DNA to RNA
                       (DNAAlphabet{4}, RNAAlphabet{2}),
                       (RNAAlphabet{4}, DNAAlphabet{2}), # RNA to DNA
                       (RNAAlphabet{2}, DNAAlphabet{4})]

    @eval function (::Type{LongSequence{$alpha}})(seq::LongSequence{$alphb})
        newseq = LongSequence{$alpha}(length(seq))
        for (i, x) in enumerate(seq)
            unsafe_setindex!(newseq, x, i)
        end
        return newseq
    end
end

for (alpha, alphb) in [(DNAAlphabet{2}, RNAAlphabet{2}),
                       (RNAAlphabet{2}, DNAAlphabet{2}),
                       (DNAAlphabet{4}, RNAAlphabet{4}),
                       (RNAAlphabet{4}, DNAAlphabet{4})]

    @eval function (::Type{LongSequence{$alpha}})(seq::LongSequence{$alphb})
        return  LongSequence{$alpha}(copy(seq.data), length(seq))
    end
end


# Concatenate multiple sequences
function LongSequence{A}(chunks::LongSequence{A}...) where {A}
    len = 0
    for chunk in chunks
        len += length(chunk)
    end
    seq = LongSequence{A}(len)
    offset = 1
    for chunk in chunks
        copyto!(seq, offset, chunk, 1, length(chunk))
        offset += length(chunk)
    end
    return seq
end

function Base.repeat(chunk::LongSequence{A}, n::Integer) where {A}
    seq = LongSequence{A}(length(chunk) * n)
    offset = 1
    for i in 1:n
        copyto!(seq, offset, chunk, 1, length(chunk))
        offset += length(chunk)
    end
    return seq
end

# Concatenation and Base.repeat operators
Base.:*(chunk::LongSequence{A}, chunks::LongSequence{A}...) where {A} =
    LongSequence{A}(chunk, chunks...)
Base.:^(chunk::LongSequence, n::Integer) = repeat(chunk, n)

function Base.similar(seq::LongSequence{A}, len::Integer = length(seq)) where {A}
    return LongSequence{A}(len)
end
