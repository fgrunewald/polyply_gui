import os
import io
import contextlib
from polyply import gen_itp
from pathlib import Path

def run_gen_itp(graph_path, outpath, force_field):
    """
    Ugly workaround calss to call gen_itp from
    the polyply library, capture the standard output
    and return it.
    """

    class input_polyply:

        def __init__(self, **kwargs):
            for key, value in kwargs.items():
                setattr(self, key, value)

    args = input_polyply(name="polyply-gui",
                         inpath=None,
                         verbosity=0,
                         seq=None,
                         seq_file=Path(graph_path),
                         outpath=Path(outpath),
                         lib=[force_field])

    with contextlib.redirect_stderr(io.StringIO()) as output:
        gen_itp(args)

    return output.getvalue().split('\n')
