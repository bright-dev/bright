"""C++ wrapper for isoname library."""
from libcpp.map cimport map
from libcpp.set cimport set

cimport std

cdef extern from "isoname.h" namespace "isoname":
    # Conversion dictionaries
    map[std.string, int] LLzz
    map[int, std.string] zzLL

    # Elemental string sets
    set[std.string] LAN
    set[std.string] ACT
    set[std.string] TRU
    set[std.string] MA
    set[std.string] FP

    # Elemental integer sets
    set[int] lan
    set[int] act
    set[int] tru
    set[int] ma
    set[int] fp

    # Current Form
    std.string CurrentForm(int)
    std.string CurrentForm(std.string)

    # LLAAAM_2_* Functions
    int LLAAAM_2_zzaaam(std.string)
    int LLAAAM_2_MCNP(std.string)

    # zzaaam_2_* Functions
    std.string zzaaam_2_LLAAAM(int)
    int zzaaam_2_MCNP(int)

    # MCNP_2_* Functions 
    int MCNP_2_zzaaam(int)
    std.string MCNP_2_LLAAAM(int)

    # mixed_2_*_ Functions
    int mixed_2_zzaaam(std.string)
    int mixed_2_zzaaam(int)
    std.string mixed_2_LLAAAM(std.string)
    std.string mixed_2_LLAAAM(int)
    int mixed_2_MCNP(std.string)
    int mixed_2_MCNP(int)

    # Helper Functions
    double nuc_weight_zzaaam(int)
    double nuc_weight(int)
    double nuc_weight(std.string)

