"""
    SobolSample(R::RandomizationMethod = NoRand()) <: DeterministicSamplingAlgorithm

Samples taken from Sobol's base-2 sequence.
"""
Base.@kwdef @concrete struct SobolSample <: DeterministicSamplingAlgorithm
    R::RandomizationMethod = NoRand()
end

function sample(n::Integer, d::Integer, S::SobolSample, T = Float64)
    @assert isinteger(log(2, n)) "Sobol nets exist only for sample sizes that are powers of 2."
    s = Sobol.SobolSeq(zeros(T, d), ones(T, d))
    return randomize(reduce(hcat, [next!(s) for i in 1:n]), S.R)
end
