__all__ = ['pmap', 'pmap_reduce']

"""Ref: Kahler: An Implementation of Discrete Exterior Calculus on Hermitian Manifolds
    code: https://github.com/aeftimia/kahler
"""
from  multiprocessing import cpu_count, Queue, Process
from operator import add

def spawn_reduce(f, fun):
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
    def func(q_in,q_out):
        while True:
            i,x = q_in.get()
            if i is None:
                break
            q_out.put((i,f(x)))
    return func

def pmap(f, X, nprocs = cpu_count()):
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