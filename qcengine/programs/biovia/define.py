import itertools as it
from typing import Any, Dict, Optional, Tuple

from qcelemental.models import AtomicInput
from qcengine.exceptions import InputError
from qcengine.util import execute

from .methods import KEYWORDS, METHODS


def define(input_model: AtomicInput) -> Tuple[bool, Dict[str, Any], Dict[str, str]]:
    """Run define and return control"""
    molecule = input_model.molecule
    charge = molecule.molecular_charge
    multiplicity = molecule.molecular_multiplicity
    method = input_model.model.method
    basis = input_model.model.basis
    keywords = input_model.keywords

    coord_str, _ = input_model.molecule.to_string(dtype="turbomole", return_data=True)

    define_input, _ = prepare_stdin(method, basis, keywords, charge, multiplicity)
    infiles = {'define.in': define_input, 'coord': coord_str}
    outfiles = ['control', 'mos', 'basis']

    commands = "ml medsci turbomole && define < define.in"

    exe_success, proc = execute(
        command=commands.split(),
        infiles=infiles,
        outfiles=outfiles,
    )

    infiles.update(proc['outfiles'])

    return exe_success, proc, infiles


def execute_define(
    inputs: Dict[str, Any],
    scratch_name: Optional[str] = None,
    timeout: Optional[int] = None,
) -> str:
    """Call define with the input define."""
    commands = "ml medsci turbomole && define < define.in".split()
    infiles = inputs["infiles"]
    outfiles = ['control']

    exe_success, proc = execute(
        command=commands.split(),
        infiles=infiles,
        outfiles=outfiles,
    )


def prepare_stdin(
    method: str, basis: str, keywords: Dict[str, Any], charge: int, mult: int, geoopt: Optional[str] = ""
) -> str:
    """Prepares a str that can be sent to define to produce the desired
    input for Turbomole."""

    # Load data from keywords
    unrestricted = keywords.get("unrestricted", False)
    grid = keywords.get("grid", "m3")

    methods_flat = list(it.chain(*[m for m in METHODS.values()]))
    if method not in methods_flat:
        raise InputError(f"Method {method} not in supported methods " f"{methods_flat}!")

    # This variable may contain substitutions that will be made to
    # the control file after it was created from a define call, e.g.
    # setting XC functionals that aren't hardcoded in define etc.
    subs = None

    def occ_num_mo_data(charge: int, mult: int, unrestricted: Optional[bool] = False) -> str:
        """Handles the 'Occupation Number & Molecular Orbital' section
        of define. Sets appropriate charge and multiplicity in the
        system and decided between restricted and unrestricted calculation.

        RHF and UHF are supported. ROHF could be implemented later on
        by using the 's' command to list the available MOs and then
        close the appropriate number of MOs to doubly occupied MOs
        by 'c' by comparing the number of total MOs and the desired
        multiplicity."""

        # Do unrestricted calculation if explicitly requested or mandatory
        unrestricted = unrestricted or (mult != 1)
        unpaired = mult - 1
        charge = int(charge)

        occ_num_mo_data_stdin = f"""eht
        y
        {charge}
        y
        """
        if unrestricted:
            # Somehow Turbomole/define asks us if we want to write
            # natural orbitals... we don't want to.
            occ_num_mo_data_stdin = f"""eht
            y
            {charge}
            n
            u {unpaired}
            *
            n
            """
        return occ_num_mo_data_stdin

    def set_method(method, grid):
        if method == "hf":
            method_stdin = ""
        elif method in METHODS["ricc2"]:
            # Setting geoopt in $ricc2 will make the ricc2 module to produce
            # a gradient.
            # Drop the 'ri'-prefix of the method string.
            geoopt_stdin = f"geoopt {method[2:]} ({geoopt})" if geoopt else ""
            method_stdin = f"""cc
                               freeze
                               *
                               cbas
                               *
                               ricc2
                               {method}
                               list models

                               {geoopt_stdin}
                               list geoopt

                               *
                               *
                            """
        elif method in METHODS["dft_hardcoded"]:
            method_stdin = f"""dft
                               on
                               func
                               {method}
                               grid
                               {grid}

                            """
        # TODO: Handle xcfuncs that aren't defined in define, e.g.
        # new functionals introduced in 7.4 from libxc. ...
        # Maybe the best idea would be to not set the functional here
        # but just turn on DFT and add it to the control file later on.
        elif method in METHODS["dft_libxc"]:
            raise InputError("libxc functionals are not supported right now.")
        return method_stdin

    # Resolution of identity
    def set_ri(keywords):
        # TODO: senex/RIJCOSX?
        ri_kws = {ri_kw: keywords.get(ri_kw, False) for ri_kw in KEYWORDS["ri"]}
        ri_stdins = {"rijk": "rijk\non\n\n", "ri": "ri\non\n\n", "marij": "marij\n\n"}
        ri_stdin = "\n".join([ri_stdins[ri_kw] for ri_kw, use in ri_kws.items() if use])
        return ri_stdin

        # ri_stdin = ""
        # # Use either RIJK or RIJ if requested.
        # if ri_kws["rijk"]:
        # ri_stdin = """rijk
        # on

        # """
        # elif ri_kws["rij"]:
        # ri_stdin = """rij
        # on

        # """
        # # MARIJ can be used additionally.
        # if ri_kws["marij"]:
        # ri_stdin += """marij

        # """
        # return ri_stdin

    # Dispersion correction
    def set_dsp(keywords):
        # TODO: set_ri and set_dsp are basically the same funtion. Maybe
        # we could abstract this somehow?
        dsp_kws = {dsp_kw: keywords.get(dsp_kw, False) for dsp_kw in KEYWORDS["dsp"]}
        dsp_stdins = {"d3": "dsp\non\n\n", "d3bj": "dsp\nbj\n\n"}
        dsp_stdin = "\n".join([dsp_stdins[dsp_kw] for dsp_kw, use in dsp_kws.items() if use])
        return dsp_stdin

    kwargs = {
        "init_guess": occ_num_mo_data(charge, mult, unrestricted),
        "set_method": set_method(method, grid),
        "ri": set_ri(keywords),
        "dsp": set_dsp(keywords),
        "title": "QCEngine Turbomole",
        "scf_conv": 8,
        "scf_iters": 150,
        "basis": basis,
    }

    stdin = """
    {title}
    a coord
    *
    no
    b
    all {basis}
    *
    {init_guess}
    {set_method}
    {ri}
    {dsp}
    scf
    conv
    {scf_conv}
    iter
    {scf_iters}

    *
    """.format(
        **kwargs
    )

    return stdin, subs
