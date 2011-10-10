"""Python header for enrichment library."""
cimport fccomp
cimport cpp_enrichment

cdef class EnrichmentParameters:
    cdef cpp_enrichment.EnrichmentParameters * ep_pointer


cdef class Enrichment(fccomp.FCComp):
    cdef cpp_enrichment.Enrichment * e_pointer


