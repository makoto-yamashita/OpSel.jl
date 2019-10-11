using OpSel
using Test

@testset "OpSel.jl" begin
    # Write your own tests here.
    OpSel.testUnequal(2045);
    OpSel.testEqual(200,50);
end
