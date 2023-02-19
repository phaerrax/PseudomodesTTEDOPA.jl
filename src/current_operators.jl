# Spin current operators 
# ======================

function J⁺tag(::SiteType"S=1/2", left_site::Int, i::Int)
  # TODO: still useful?
  if i == left_site
    str = "σx"
  elseif i == left_site + 1
    str = "σy"
  else
    str = "Id"
  end
  return str
end
function J⁺tag(::SiteType"vecS=1/2", left_site::Int, i::Int)
  if i == left_site
    str = "vecσx"
  elseif i == left_site + 1
    str = "vecσy"
  else
    str = "vecId"
  end
  return str
end
function J⁺tag(::SiteType"HvS=1/2", left_site::Int, i::Int)
  return J⁺tag(SiteType("vecS=1/2"), left_site, i)
end

function J⁻tag(::SiteType"S=1/2", left_site::Int, i::Int)
  # Just as `J⁺tag`, but for σʸ⊗σˣ
  if i == left_site
    str = "σy"
  elseif i == left_site + 1
    str = "σx"
  else
    str = "Id"
  end
  return str
end
function J⁻tag(::SiteType"vecS=1/2", left_site::Int, i::Int)
  if i == left_site
    str = "vecσy"
  elseif i == left_site + 1
    str = "vecσx"
  else
    str = "vecId"
  end
  return str
end
function J⁻tag(::SiteType"HvS=1/2", left_site::Int, i::Int)
  return J⁻tag(SiteType("vecS=1/2"), left_site, i)
end

function spin_current_op_list(sites::Vector{Index{Int64}})
  N = length(sites)
  # Check if all sites are spin-½ sites.
  if all(x -> SiteType("S=1/2") in x, sitetypes.(sites))
    st = SiteType("S=1/2")
    MPtype = MPO
  elseif all(x -> SiteType("vecS=1/2") in x, sitetypes.(sites))
    st = SiteType("vecS=1/2")
    MPtype = MPS
  elseif all(x -> SiteType("HvS=1/2") in x, sitetypes.(sites))
    st = SiteType("HvS=1/2")
    MPtype = MPS
  else
    throw(ArgumentError("spin_current_op_list works with SiteTypes "*
                        "\"S=1/2\", \"vecS=1/2\" or \"HvS=1/2\"."))
  end
  #
  J⁺ = [MPtype(sites, [J⁺tag(st, k, i) for i = 1:N]) for k = 1:N-1]
  J⁻ = [MPtype(sites, [J⁻tag(st, k, i) for i = 1:N]) for k = 1:N-1]
  return -0.5 .* (J⁺ .- J⁻)
end
