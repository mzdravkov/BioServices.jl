module KEGG

import HTTP

const base_url = "https://rest.kegg.jp/"


@enum KEGG_DB begin
  PATHWAY
  BRITE
  MODULE
  ORTHOLOGY
  GENES
  VIRAL_GENES
  VIRAL_PEPTIDES
  ADDENDUM_GENES
  GENOME
  COMPOUND
  GLYCAN
  REACTION
  RCLASS
  ENZYME
  NETWORK
  VARIANT
  DISEASE
  DRUG
  DGROUP
  ORGANISM
end


struct PathwayHeader
  id::String
  name::String
end

struct BriteHeader
  id::String
  name::String
end

struct ModuleHeader
  id::String
  name::String
end

struct KOHeader
  id::String
  symbols::AbstractArray{String}
  name::String
  enzymes::AbstractArray{String}
end

struct GenomicLocation
  chr::String
  strand::Union{String, Nothing}
  start_pos::Union{Integer, Nothing}
  end_pos::Union{Integer, Nothing}
end

struct GeneHeader
  id::String
  type::String
  locations::AbstractArray{GenomicLocation}
  symbols::AbstractArray{String}
  description::String
end

struct ViralGeneHeader
  id::String
  symbols::AbstractArray{String}
  organism::String
  description::String
end

struct ViralPeptideHeader
  id::String
  description::String
end

struct AddendumGeneHeader
  id::String
  symbol::Union{Nothing, String}
  name::String
  enzymes::AbstractArray{String}
end


function parse_ko_header(line)
  id, rest = split(line, "\t")
  if !occursin(";", rest)
    return KOHeader(id, [], rest, [])
  end
  symbols, rest = split(rest, "; ")
  symbols = split(symbols, ", ")
  if occursin(" [EC:", rest)
    name, rest = split(rest, " [EC:")
    enzymes, _ = split(rest, "]")
    enzymes = split(enzymes, " ")
    KOHeader(id, symbols, name, enzymes)
  else
    KOHeader(id, symbols, rest, [])
  end
end


function parse_gene_header(line)
  id, type, location_col, gene_info = split(line, "\t")

  # The gene_info column would usually have the form
  # symbol1, symbolN; desciption
  # However, some genes have no HGNC symbols, so we
  # assume that in abscence of a semicolon, the whole
  # column is a description.
  if occursin("; ", gene_info)
    symbols, description = split(gene_info, "; ")
    symbols = split(symbols, ", ")
  else
    symbols = []
    description = gene_info
  end

  locations = []
  for location_str in split(location_col, ", ")
    # For some genes, there's only information about
    # the chromosome, but not a precise position on it.
    if !occursin(":", location_str)
      push!(locations, GenomicLocation(location_str,
                                       nothing,
                                       nothing,
                                       nothing))
      continue
    end

    chr, chr_location = split(location_str, ":")
    if startswith(chr_location, "complement")
      strand = "+"
      m = match(r"complement\(([0-9\.]+)\)", chr_location)[1]
      start_pos, end_pos = split(m, "..")
    else
      strand = "-"
      start_pos, end_pos = split(chr_location, "..")
    end
    location = GenomicLocation(chr,
                               strand,
                               parse(Int64, start_pos),
                               parse(Int64, end_pos))
    push!(locations, location)
  end

  GeneHeader(id, type, locations, symbols, nothing, description)
end

function parse_viral_gene_header(line)
  id, gene_info = split(line, "\t")

  # The gene_info column would usually have the form
  # symbol1, symbolN; organism; desciption
  symbols, organism, description = split(gene_info, "; ")
  symbols = split(symbols, ", ")

  ViralGeneHeader(id, symbols, organism, description)
end


function parse_addendum_gene_header(line)
  id, gene_info = split(line, "\t")

  # The gene_info column would usually have the form
  # symbol?; name <optional enzymes list>
  if occursin("; ", gene_info)
    symbol, rest = split(gene_info, "; ")
  else
    symbol = nothing
    rest = gene_info
  end

  if occursin(" (EC:", rest)
    name, rest = split(rest, " (EC:")
    enzymes, _ = split(rest, ")")
    enzymes = split(enzymes, " ")
    AddendumGeneHeader(id, symbol, name, enzymes)
  else
    AddendumGeneHeader(id, symbol, rest, [])
  end
end


function parse_list_line(line, database)
  if database == PATHWAY
    return PathwayHeader(split(line, "\t")...)
  elseif database == BRITE
    return BriteHeader(split(line, "\t")...)
  elseif database == MODULE
    return ModuleHeader(split(line, "\t")...)
  elseif database == ORTHOLOGY
    return parse_ko_header(line)
  elseif database == GENES
    return parse_gene_header(line)
  elseif database == VIRAL_GENES
    return parse_viral_gene_header(line)
  elseif database == VIRAL_PEPTIDES
    return ViralPeptideHeader(split(line, "\t")...)
  elseif database == ADDENDUM_GENES
    return parse_addendum_gene_header(line)
  else
    throw(MethodError(parse_list_line, "Parsing list response not implemented for database \"$database\""))
  end
end


function parse_list(response, database)
  lines = eachline(IOBuffer(response.body))
  [parse_list_line(line, database) for line in lines]
end


function list_endpoint(url, db)
  res = HTTP.request("GET", url)
  return parse_list(res, db)
end


function list_pathways(organism::Union{Nothing, String}=nothing)
  url = string(base_url, "list/pathway")
  if organism !== nothing
    url = string(url, "/", organism)
  end

  return list_endpoint(url, PATHWAY)
end

function list_brite(option::Union{Nothing, String}=nothing)
  url = string(base_url, "list/brite")
  if option !== nothing
    url = string(url, "/", option)
  end

  return list_endpoint(url, BRITE)
end

function list_modules()
  url = string(base_url, "list/module")
  return list_endpoint(url, MODULE)
end


function list_orthologs()
  url = string(base_url, "list/orthology")
  return list_endpoint(url, ORTHOLOGY)
end


"""
list_genes returns a list of genes from a particular database.

The possible gene databases in KEGG are:
- <org>: Genes for the specified organism.
- vg: Genes in viruses category
- vp: Mature peptides in viruses
- ag: Genes in addendum category
"""
function list_genes(db::String)
  url = string(base_url, "list/", db)
  kegg_db = if db == "vg"
    VIRAL_GENES
  elseif db == "vp"
    VIRAL_PEPTIDES
  elseif db == "ag"
    ADDENDUM_GENES
  else
    GENES
  end
  list_endpoint(url, kegg_db)
end


function list_genomes()
  url = string(base_url, "list/genome")
  return list_endpoint(url, GENOME)
end


function list_compounds()
  url = string(base_url, "list/compound")
  return list_endpoint(url, COMPOUND)
end


function list_glycans()
  url = string(base_url, "list/glycan")
  return list_endpoint(url, GLYCAN)
end

# function kegg_list(database::Union{KEGG_DB, String};
#                    organism::Union{Nothing, String}=nothing,
#                    option::Union{Nothing, String}=nothing)
#   if organism != nothing && database ∉ [PATHWAY BRITE]
#     throw(MethodError(kegg_list, "The organism parameter is only available for the \"pathway\" and \"brite\" databases."))
#   end

#   if isa(database, KEGG_DB)
#     database_name = string(database) |> lowercase
#   else
#     database_name = database
#   end

#   url = string(base_url, "list/", database_name)
#   if organism !== nothing
#     url = string(url, "/", organism)
#   elseif option !== nothing
#     url = string(url, "/", option)
#   end

#   res = HTTP.request("GET", url)
#   return parse_list(res, database)
# end

end
