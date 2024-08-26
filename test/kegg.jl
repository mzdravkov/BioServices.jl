@testset "KEGG" begin
  @testset "list" begin
    pathways = KEGG.kegg_list(KEGG.PATHWAY)
    println(pathways)
  end
end
