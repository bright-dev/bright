template = """\
[CHAR] LWR Fuel Cell
C 
C ********************
C **** Cell Cards ****
C ********************
C 
1  1 -{FuelDensity}   -1 -3    imp:n=1 vol={FuelCellVolume}  $Fuel Region
2  2 -{CladDensity}    1 -2 -3 imp:n=1           $Clad Region
3  3 -{CoolDensity}    2 -3    imp:n=1           $Water Region
10 0  3               imp:n=0           $Void Region

C 
C ***********************
C **** Surface Cards ****
C ***********************
C 
1 CZ {FuelCellRadius} $Fuel Cylinder
2 CZ {CladCellRadius} $Clad Cylinder
*3 BOX 0 0 0 {UnitCellHalfPitch} 0 0 0 {UnitCellHalfPitch} 0 0 0 {UnitCellHeight} $box of interest

C 
C *******************
C **** Burn Card ****
C *******************
C 
{Burn}
C 
C ************************
C **** Material Cards ****
C ************************
C 
{Mat1}
m2 40090 1.0                   $Cladding
m3 01001 0.66667 08016 0.33333 $Light water
mt3 lwtr.04t                   $Scattering at 600 K
{MatOther}
C 
C ****************************
C **** Perturbation Cards ****
C ****************************
C 
{Pert}
C 
C ***********************
C **** Physics Cards ****
C ***********************
C 
c phys:n 10.
c mode n
C
C **********************
C **** Source Cards ****
C **********************
C
c kcode {kParticles} 1.0 {kCyclesSkip} {kCycles} 
c ksrc 0.1 0.1 5
c
SDEF POS=0.1 0.1 5
NPS {kParticles}
SI1 {GroupStructure}
SP1 0 0.2 0.2 0.2 0.2 0.2
c imp:n 1 1 1 0
C 
C *********************
C **** Tally Cards ****
C *********************
C
{Tally}
{TallyScat}
C 
C **********************
C **** Output Cards ****
C **********************
C 
print -10 -20 -30 -35 -40 -50 -60 -62 -72 -85 -90 -98
      -100 -102 -110 -120 -126 -128 -130 -140 -150 -160 
      -161 -162 -170 -175 -178 -180 -190 -198 -200 """


def template_extras(moffset=3):
    """Computes extra fields for the template that must be computed dymaically."""
    from char import defchar

    mat_number = {}

    mnum = moffset
    for nuc in defchar.core_load["MCNP"]:
        mnum = mnum + 1
        mat_number[nuc] = mnum

    number_mat = {v: k for k, v in mat_number.items()}

    extra_fill = {'mat_number': mat_number, 'number_mat': number_mat}

    return extra_fill
