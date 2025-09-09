__all__ = ['pmap', 'pmap_reduce']

r"""
Parallel utilities for simple multiprocessing map/reduce patterns.

This module provides lightweight wrappers around Python's multiprocessing
module to parallelize independent computations over iterables. It is used
internally to speed up tasks like building isotypic bases where many small
independent jobs can be executed concurrently.

Environment versions
- SageMath 10.4
- Python 3.11
- NumPy 1.26
- SciPy 1.14

Notes
- These functions use processes (not threads). Functions f must be picklable.
- On Windows, protect entry points with "if __name__ == '__main__':" before
  calling pmap/pmap_reduce to avoid recursive process spawning.
- The API is pure Python and does not depend on Sage; examples can be run in
  a standard Python interpreter or inside a Sage session.

REFERENCES:
- Inspired by code patterns in: Kahler — An Implementation of Discrete Exterior
  Calculus on Hermitian Manifolds (https://github.com/aeftimia/kahler)
"""
from  multiprocessing import cpu_count, Queue, Process
from operator import add

def spawn_reduce(f, fun):
    r"""
    Internal helper: create a worker function that reduces results.

    INPUT:
    - f -- function to apply to each item
    - fun -- associative binary function to combine partial results

    OUTPUT:
    - callable expecting (q_in, q_out); used as Process target
    """
    def func(q_in,q_out):
        i,x = q_in.get()
        if i is None:
            q_out.put(None)
            return
        ret = f(x)
        while True:
            i,x = q_in.get()
            if i is None:
                break
            ret = fun(ret, f(x))
        q_out.put(ret)
    return func

def spawn(f):
    r"""
    Internal helper: create a worker function for mapping.

    INPUT:
    - f -- function to apply to each item

    OUTPUT:
    - callable expecting (q_in, q_out); used as Process target
    """
    def func(q_in,q_out):
        while True:
            i,x = q_in.get()
            if i is None:
                break
            q_out.put((i,f(x)))
    return func

def pmap(f, X, nprocs = cpu_count()):
    r"""
    Parallel map over iterable X using up to nprocs worker processes.

    INPUT:
    - f -- callable applied to each element of X (must be picklable)
    - X -- iterable of inputs
    - nprocs -- number of processes (default: cpu_count())

    OUTPUT:
    - list with f(x) for x in X, preserving the original order

    EXAMPLES::

        sage: from pysymmetry.parallel import pmap
        sage: pmap(lambda t: t*t, [0,1,2,3], nprocs=2)
        [0, 1, 4, 9]
    """
    q_in   = Queue(True)
    q_out  = Queue()

    proc = [Process(target=spawn(f),args=(q_in, q_out)) for _ in range(nprocs)]
    for p in proc:
        p.daemon = True
        p.start()
        proc

    [q_in.put((i,x)) for i,x in enumerate(X)]
    [q_in.put((None,None)) for _ in range(nprocs)]

    result = [q_out.get() for _ in range(len(X))]

    [p.join() for p in proc]
    result = [x for i,x in sorted(result)]
    return result

def pmap_reduce(f, X, fun=add, nprocs = cpu_count()):
    r"""
    Parallel map followed by reduction across worker-local partial results.

    INPUT:
    - f -- callable applied to each element of X (must be picklable)
    - X -- iterable of inputs
    - fun -- associative binary function used to combine partial results
             (default: operator.add)
    - nprocs -- number of processes (default: cpu_count())

    OUTPUT:
    - single value equal to fun(f(x0), fun(f(x1), ...)) in unspecified tree order

    EXAMPLES::

        sage: from pysymmetry.parallel import pmap_reduce
        sage: # Sum of squares of 0..4
        sage: pmap_reduce(lambda t: t*t, range(5), fun=lambda a,b: a+b, nprocs=2)
        30
    """
    q_in   = Queue(True)
    q_out  = Queue()
    proc = [Process(target=spawn_reduce(f, fun),args=(q_in, q_out)) for _ in range(nprocs)]
    for p in proc:
        p.daemon = True
        p.start()
        proc

    [q_in.put((i,x)) for i,x in enumerate(X)]
    [q_in.put((None,None)) for _ in range(nprocs)]

    result = None
    for _ in range(nprocs):
        partial = q_out.get()
        if partial is not None:
            if result is None:
                result = partial
            else:
                result = fun(result, partial)

    [p.join() for p in proc]
    return result