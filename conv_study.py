from fictdom_dgT import SolveFdT
from ngsolve import TaskManager, SetNumThreads
from math import log

import argparse

SetNumThreads(8)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Solve unfitted/ficticious domain Poisson problem with Nitsche boundary value imposition and Trefftz.')
    parser.add_argument('-ex', '--example', type=int, default=1, help='Example to run, deafult=1')
    parser.add_argument('-m', '--method', default="trefftz", help='Method used: options trefftz, emb_trefftz, trefftz_mixed, l2, default: trefftz')
    parser.add_argument('-s', '--stabil', default="global", help='Patchwise or global ghost-penalties: options patch, default: global')
    parser.add_argument('-o', '--order', type=int, default=3, help='order of discretisation & mesh deformation, default: 3')
    parser.add_argument('-nr', '--n_ref', type=int, default=5, help='number of refinements, default: 5')
    parser.add_argument('-qm', '--quad_mesh', type=int, default=0, help='Use a quad mesh, enter either int 0 for false or int 1 for true. default: 0')
    parser.add_argument('-sm', '--struc_mesh', type=int, default=0, help='Use structured mesh, enter either int 0 for false or int 1 for true. default: 0')
    parser.add_argument('-def', '--deformation', type=int, default=1, help='Use ioparametric mapping for higher-order geometry approximation. default: 1')

    args = parser.parse_args()
    options = vars(args)


def conv_study(options):
    outfile_name = f'out/conv_ex{options["example"]}_{options["method"]}_gp{options["stabil"]}_p{options["order"]}_nref{options["n_ref"]}_qm{options["quad_mesh"]}_sm{options["struc_mesh"]}_def{options["deformation"]}.dat'
    f = open(outfile_name, 'w')
    f.close()

    l2errors = []
    for i in range(options['n_ref']):
        with TaskManager():
            l2err, nd, times = SolveFdT(
                example=options['example'], method=options['method'],
                maxh=0.5**(i + 1), order=options['order'],
                quad_mesh=bool(options['quad_mesh']),
                struc_mesh=bool(options['struc_mesh']),
                deformation=bool(options['deformation']))
        f = open(outfile_name, 'a')
        f.write(f'{i}\t{l2err}\t{nd}\n')
        f.close
        print(f'Single run result: {l2err} {nd}')
        l2errors.append(l2err)

    print(f'final l2 error stat: {l2errors}')
    eocs = [log(l2errors[j - 1] / l2errors[j]) / log(2)
            for j in range(1, len(l2errors))]
    print(f'eocs  : {eocs}')


if __name__ == '__main__':
    conv_study(options)
