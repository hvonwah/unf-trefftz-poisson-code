from fictdom_dgT import SolveFdT
from ngsolve import TaskManager, SetNumThreads
from math import log

import argparse

parser = argparse.ArgumentParser(description='Solve unfitted/ficticious domain Poisson problem with Nitsche boundary value imposition and Trefftz.')
parser.add_argument('-th', '--threads', type=int, default=1, help='number of shared memeory paralel threads, default: 1')
parser.add_argument('-ex', '--example', type=int, default=1, help='Example to run, deafult=1')
parser.add_argument('-s', '--stabil', default="global", help='Patchwise or global ghost-penalties: options patch, default: global')
parser.add_argument('-o', '--order', type=int, default=3, help='order of discretisation & mesh deformation, default: 3')
parser.add_argument('-hmax', '--mesh_size', type=float, default=0.5, help='Mesh size of the background mesh, default: 0.5')
parser.add_argument('-qm', '--quad_mesh', type=int, default=0, help='Use a quad mesh, enter either int 0 for false or int 1 for true. default: 0')
parser.add_argument('-sm', '--struc_mesh', type=int, default=0, help='Use structured mesh, enter either int 0 for false or int 1 for true. default: 0')
parser.add_argument('-def', '--deformation', type=int, default=1, help='Use ioparametric mapping for higher-order geometry approximation. default: 1')

args = parser.parse_args()
options = vars(args)

SetNumThreads(options['threads'])

outfile_name = f'out/timings_threads{options["threads"]}_example{options["example"]}_gp{options["stabil"]}_p{options["order"]}_h{options["mesh_size"]}_qm{options["quad_mesh"]}_sm{options["struc_mesh"]}_def{options["deformation"]}.dat'

with open(outfile_name, 'w') as fid:
    fid.write('method assemble embedding solve')

    for method in ['l2', 'trefftz']:
        with TaskManager():
            l2err, nd, times = SolveFdT(
                example=options['example'], method=method,
                maxh=options['mesh_size'], order=options['order'],
                quad_mesh=bool(options['quad_mesh']),
                struc_mesh=bool(options["struc_mesh"]),
                deformation=bool(options['deformation']))

        fid.write(f'\n{method} ')
        for key in ['assemble', 'embedding', 'solve']:
            fid.write(f'{times[key]} ')
