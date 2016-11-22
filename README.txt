

So far, the metabolite a in WS model represented both energy and processed metabolites. We seek to separate metabolism into two distinct paths: 1 - Catabolism, i.e. Cutting the sugars in order to extract ATP (From glycolysis) & RedOx potential (to eventually produce ATP from respiration). 2 - Anabolism, i.e Molecules are processed into precursor of the cell mass: amino acids, dNTPs, rNTPs, etc.

Respiration requires maintenance of the pH gradient between the exterior/interior of the cell. Exterior would be the periplasmic space for gram-positive bacteria. This is costly as it requires maintenance of the membrane potential.


Edit1: The reaction of ribosome binding to mRNA now writes:
10 * r + m_x <-> c_x @ k_b / 10, k_u

The total amount of Ribosomes R_Tot now is:
R_Tot = r + 10 \cdot \sum_x c_x + zm_x

The translation rate per ribosome is now the translation rate per 10 ribosomes and writes:
\gamma(a) = \gamma_{max} \cdot a / (K_{\gamma} + a)
with \gamma_{max} equals ten times what it was previously, i.e. now \gamma_{max} = 12600




