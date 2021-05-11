"""
General-purpose miscellaneous declarations useful to all submodules.
"""

import numpy as np
from functools import wraps

Real = (int, np.int32, np.int64, float, np.float32, np.float64)
NPReal = (*Real, np.ndarray)
Int = (int, np.int32, np.int64)
NPInt = (*Int, np.ndarray)
List_like = (list, tuple)
NPList_like = (*List_like, np.ndarray)
Iterable = (*List_like, set)
NPIterable = (*Iterable, np.ndarray)
Parameter = (*Real, type(None))
Function = type(lambda:0)

class BadType(Exception):
    def __init__(self, argname, given=None, expected=None):
        if given != None and expected != None:
            self.message = f"arg `{argname}` is not of expected type (expected {expected}, but got {given})"
        elif expected == None:
            self.message = f"arg `{argname}` is not of expected type (expected {expected})"
        elif given == None:
            self.message = f"arg `{argname}` is not of expected type (got {given})"
        else:
            self.message = f"arg `{argname}` is not of expected type"
        Exception.__init__(self, self.message)

do_typecheck = True
def typecheck(types, conds=dict()):
    """
    Decorator which performs a high-level check that all arguments passed to a
    given function it satisfy specified type requirements and additional
    conditions expressable as single-argument functions which return `True` or
    `False`.

    Parameters
    ----------
    types : dict-like[str : type or Tuple of types or None]
        Dict-like mapping argument names (as strs) to types or Tuples of types
        or `None`. If an argument name is mapped to `None`, then all types are
        accepted, and additional conditions on the argument are still performed.
    conds : dict-like[int : True/False Callable or Iterable of True/False Callables]
        Dict-like mapping argument names (as strs) to single-argument Callables
        which return `True` or `False` (with `True` being returned when a
        favorable condition is satisfied) or an Iterable of such Callables.
    """
    def decorator(f):
        @wraps(f)
        def checker(*args, **kwargs):
            if do_typecheck:
                argnames = f.__code__.co_varnames
                argvs = dict()
                argvs.update(kwargs)
                argvs.update(dict(zip(argnames, args)))

                for arg in argvs.keys():
                    typespec = types.get(arg, None)
                    if typespec == None:
                        pass
                    elif not isinstance(argvs[arg], typespec):
                        if isinstance(typespec, tuple):
                            expected = [t.__name__ for t in typespec]
                        else:
                            expected = typespec.__name__
                        raise BadType(arg, type(argvs[arg]).__name__, expected)

                    conditions = conds.get(arg, list())
                    if isinstance(conditions, type(lambda:None)):
                        conditions = [conditions]
                    elif not isinstance(conditions, (tuple, set, list)):
                        raise Exception(f"{f.__name__}: typecheck: invalid conditions provided for arg `{arg}`")
                    for cond in conditions:
                        try:
                            cond_check = cond(argvs[arg])
                        except Exception as e:
                            raise Exception(f"{f.__name__}: typecheck: error occurred while testing conditions for arg {arg}:\n{e}")
                        if not cond_check:
                            raise Exception(f"{f.__name__}: typecheck: arg `{arg}` did not meet specified conditions")
            return f(*args, **kwargs)
        return checker
    return decorator

def gen_table_fmt(label_fmts, s="  ", L=12, P=5, K=2) -> (str, str):
    """
    Generate the column labels and format string of a table from a list of
    tuples following
        (
            'column label',
            x in {'s','s>','i','f','g','e'},
            {l: length override, p: precision override} (optional)
        )
    """
    head = ""
    lines = ""
    fmt = ""
    names = list()
    for label_fmt in label_fmts:
        names.append(label_fmt[0])
        overrides = dict() if len(label_fmt) < 3 else label_fmt[2]
        l = overrides.get("l",
            max(int((len(label_fmt[0])+K-1)/K)*K, L*(label_fmt[1] in ['e','f','g']))
        )
        p = overrides.get("p",
            l-7 if (l-7 >= 1 and l-7 <= P) else P
        )
        head += "{:"+str(l)+"s}"+s
        lines += l*"-" + s
        if label_fmt[1] == 's':
            fmt += "{:"+str(l)+"s}"+s
        elif label_fmt[1] == 's>':
            fmt += "{:>"+str(l)+"s}"+s
        elif label_fmt[1] == 'i':
            fmt += "{:"+str(l)+".0f}"+s
        elif label_fmt[1] in ['e', 'f', 'g']:
            fmt += "{:"+str(l)+"."+str(p)+label_fmt[1]+"}"+s
        else:
            raise Exception("Format is not one of {'s', 's>', 'i', 'f', 'g', 'e'}")
    head = head[:-len(s)]
    lines = lines[:-len(s)]
    fmt = fmt[:-len(s)]
    return head.format(*names)+"\n"+lines, fmt

def print_write(outfile, s, end="\n", flush=True) -> None:
    print(s, end=end, flush=flush)
    outfile.write(s+end)
    if flush:
        outfile.flush()
    return None

