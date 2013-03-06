

cdef class EnrichmentParameters:
    """This class is a collection of values that mirror the attributes in 
    Enrichment that are required for the cascade model to run. Like 
    ReactorParameters, this class takes no arguments on initialization.  An 
    empty ErichmentParameters instance has all values set to zero."""

    def __cinit__(self):
        self.ep_pointer = new cpp_enrichment.EnrichmentParameters()

    def __dealloc__(self):
        del self.ep_pointer


    #
    # Class Attributes
    #

    property alpha_0:
        """The :math:`\\alpha_0` attribute specifies the overall stage separation factor
        for the cascade.  This should be set on initialization.  Values should be
        greater than one.  Values less than one represent de-enrichment."""
        def __get__(self):
            return self.ep_pointer.alpha_0

        def __set__(self, value):
            self.ep_pointer.alpha_0 = <double> value


    property Mstar_0:
        """The :math:`M^*_0` represents a first guess at what the `Mstar` should be.
        The value of Mstar_0 on initialization should be in the ballpark
        of the optimized result of the Mstar attribute.  However, :math:`M^*_0` must
        always have a value between the weights of the j and k key components."""
        def __get__(self):
            return self.ep_pointer.Mstar_0

        def __set__(self, value):
            self.ep_pointer.Mstar_0 = <double> value


    property j:
        """This is an integer in zzaaam-form that represents the jth key component.
        This nuclide is preferentially enriched in the product stream.
        For standard uranium cascades j is 922350 (ie U-235)."""
        def __get__(self):
            return self.ep_pointer.j

        def __set__(self, value):
            self.ep_pointer.j = nucname.zzaaam(value)


    property k:
        """This is an integer in zzaaam-form that represents the kth key component.
        This nuclide is preferentially enriched in the waste stream.
        For standard uranium cascades k is 922380 (ie U-238)."""
        def __get__(self):
            return self.ep_pointer.k

        def __set__(self, value):
            self.ep_pointer.k = nucname.zzaaam(value)


    property N0:
        """This is the number of enriching stages initially guessed by the user."""
        def __get__(self):
            return self.ep_pointer.N0

        def __set__(self, value):
            self.ep_pointer.N0 = <double> value


    property M0:
        """This is the number of stripping stages initially guessed by the user."""
        def __get__(self):
            return self.ep_pointer.M0

        def __set__(self, value):
            self.ep_pointer.M0 = <double> value


    property xP_j:
        """This is the target enrichment of the jth isotope in the
        product stream mat_prod.  The :math:`x^P_j` value is set by 
        the user at initialization or run-time.  For typical uranium 
        vectors, this value is about U-235 = 0.05."""
        def __get__(self):
            return self.ep_pointer.xP_j

        def __set__(self, value):
            self.ep_pointer.xP_j = <double> value


    property xW_j:
        """This is the target enrichment of the jth isotope in the
        waste stream ms_tail.  The :math:`x^W_j` value is set by the 
        user at initialization or runtime.  For typical uranium vectors,
        this value is about U-235 = 0.0025."""
        def __get__(self):
            return self.ep_pointer.xW_j

        def __set__(self, value):
            self.ep_pointer.xW_j = <double> value



desc = {
    'docstrings': {
        'module': """Python wrapper for enrichment parameters.""", 
        'class': class_ds,
        'attrs': {
            'sepeff_ds': sepeff_ds,
            },
        'methods': {
            'initialize': init_ds, 
            'calc_params': calc_params_ds,
            'calc': calc_ds,
            },
        },
    'attrs': {
        'sepeff': 'sepeff_t',
        },
    'methods': {
        ('Reprocess', ('sepeff', 'sepeff_t'), ('name', 'str', '""')): None,
        },
    }
