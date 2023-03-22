# ------------------------------ LOAD LIBRARIES -------------------------------
from netgen.geom2d import SplineGeometry
from ngsolve import *
from ngsolve.internal import *
from xfem import *
from xfem.lsetcurv import *
from ngstrefftz import *
import argparse

ngsglobals.msg_level = 2
SetNumThreads(4)

parser = argparse.ArgumentParser(description='Solve unfitted/ficticious domain convection-diffusion problem with Nitsche to impose dirchelt boundary value imposition and usimg eithe the DG or embedded Trefftz DG method.')
parser.add_argument('-m', '--method', default='emb_trefftz', help='Method used: options: emb_trefftz, l2, default: emb_trefftz')
parser.add_argument('-o', '--order', type=int, default=4, help='order of discretisation & mesh deformation, default: 4')
parser.add_argument('-nr', '--n_ref', type=int, default=2, help='number of refinements, default: 2')
parser.add_argument('-qm', '--quad_mesh', type=int, default=0, help='Use a quad mesh, enter either int 0 for false or int 1 for true. default: 0')
parser.add_argument('-sm', '--struc_mesh', type=int, default=0, help='Use structured mesh, enter either int 0 for false or int 1 for true. default: 0')
parser.add_argument('-def', '--deformation', type=int, default=1, help='Use ioparametric mapping for higher-order geometry approximation. default: 1')

args = parser.parse_args()
options = vars(args)
print(options)


# -------------------------------- PARAMETERS ---------------------------------
# Quadrilateral (or simplicial mesh)
quad_mesh = bool(options['quad_mesh'])
# Mesh diameter
maxh = 0.2 * 0.5**options['n_ref']
# Finite element space order
order = options['order']

# Geometry and Mesh
square = SplineGeometry()
square.AddRectangle((-1, -1), (1, 1), bcs=['bottom', 'right', 'top', 'left'])
ngmesh = square.GenerateMesh(maxh=maxh, quad_dominated=quad_mesh)
mesh = Mesh(ngmesh)
h = specialcf.mesh_size

# Diffusion parameter
alpha = 0.001

# Stabilization parameter for ghost-penalty
gamma_stab = 0.001 * (alpha / h**2 + 1 / h)
# Stabilization parameter for Nitsche
stab_nitsche = 10 * alpha / h * order * order
# Stabilization parameter for interior penalty
stab_dg = 10 * alpha / h * order * order


# ------------------------------- PROBLEM DATA --------------------------------
R = 0.25
r2 = x**2 + y**2
r4 = r2**2
r = sqrt(r2)
levelset = R - r

w = CF((1 + R**2 * (y**2 - x**2) / r4, -2 * R**2 * x * y / r4))

Draw(w, mesh, "w")

exact = CF(0)  # unknown
coef_f = 0
uD = 1

# ---------------- LEVELSET APPROXIMATION AND ELEMENT MARKINGS ----------------
# Higher order level set approximation
lsetmeshadap = LevelSetMeshAdaptation(mesh, order=order, threshold=0.1,
                                      discontinuous_qn=True)
deformation = lsetmeshadap.CalcDeformation(levelset)
lsetp1 = lsetmeshadap.lset_p1

# Element, facet and dof marking w.r.t. boundary approximation with lsetp1:
ci = CutInfo(mesh, lsetp1)
els_hasneg = ci.GetElementsOfType(HASNEG)
els_if = ci.GetElementsOfType(IF)
# Facets for ghost penalty stabilization
facets_gp = GetFacetsWithNeighborTypes(mesh, a=els_hasneg, b=els_if)
# Facets with parts inside the domain
facets_dg = GetFacetsWithNeighborTypes(mesh, a=els_hasneg, b=els_hasneg)


# --------------------------- FINITE ELEMENT SPACE ----------------------------
Vhbase = monomialfespace(mesh, order=order, dirichlet=[], dgjumps=True)
Vh = Restrict(Vhbase, els_hasneg)

if options['method'] == 'emb_trefftz':
    Vhbase_om2 = monomialfespace(mesh, order=order - 2, dirichlet=[],
                                 dgjumps=True)
    Vh_om2 = Restrict(Vhbase_om2, els_hasneg)

gfu = GridFunction(Vh)

# -------------------------- VARIATIONAL FORMULATION --------------------------
u, v = Vh.TnT()
nh = 1.0 / Norm(grad(lsetp1)) * grad(lsetp1)
nF = specialcf.normal(mesh.dim)
wn = w * nF
flux_u = -0.5 * (grad(u) + grad(u.Other())) * nF
flux_v = -0.5 * (grad(v) + grad(v.Other())) * nF
jump_u = u - u.Other()
jump_v = v - v.Other()


dX = dx(definedonelements=els_hasneg, deformation=deformation)
# Element-wise integrals
dx = dCut(lsetp1, NEG, definedonelements=els_hasneg, deformation=deformation)
# Interior skeleton integrals:
dk = dCut(lsetp1, NEG, skeleton=True, definedonelements=facets_dg,
          deformation=deformation)
ddk = dCut(lsetp1, NEG, element_boundary=True, definedonelements=els_hasneg,
           deformation=deformation)
db = ds(skeleton=True, deformation=deformation,
        definedon=mesh.Boundaries("left"))
# Domain boundary integrals
ds = dCut(lsetp1, IF, definedonelements=els_if, deformation=deformation)
# Ghost penalty integrals
dw = dFacetPatch(definedonelements=facets_gp, deformation=deformation)

a = RestrictedBilinearForm(Vh, "a", els_hasneg, facets_dg, check_unused=False)
a += alpha * grad(u) * grad(v) * dx
a += stab_dg * jump_u * jump_v * dk
a += alpha * (flux_u * jump_v + flux_v * jump_u) * dk
a += (stab_dg * u * v + alpha * (grad(u) * nF * v + grad(v) * nF * u)) * db
a += alpha * (-grad(u) * nh * v - grad(v) * nh * u) * ds
a += stab_nitsche * u * v * ds

a += w * grad(u) * v * dx
a += wn * IfPos(wn, 0, u.Other() - u) * v * ddk

a += gamma_stab * (u - u.Other()) * (v - v.Other()) * dw

f = LinearForm(Vh)
f += coef_f * v * dx
f += (- alpha * grad(v) * nh * uD + stab_nitsche * uD * v) * ds

# Assemble DG system
with TaskManager():
    f.Assemble()
    a.Assemble()


# ----------------------------- SET-UP EMBEDDING ------------------------------
if options['method'] == 'emb_trefftz':
    def Lap(u):
        return sum(Trace(u.Operator('hesse')))

    eps = 1e-7

    vt = Vh_om2.TestFunction()
    op = (-alpha * Lap(u) + w * grad(u)) * vt * dx
    lop = coef_f * vt * dx
    with TaskManager():
        PP, uf = TrefftzEmbedding(bf=op, fes=Vh, lf=lop, test_fes=Vh_om2)
    ncutels = sum(els_hasneg)
    print(f'n neg els: {ncutels}')

    print(f'PP.shape: {PP.shape})', end='')
    print(f'per element: {PP.shape[0] / ncutels} -> {PP.shape[1] / ncutels}')

    PPT = PP.CreateTranspose()
    TA = PPT @ a.mat @ PP


# ------------------------------- SOLVE PROBLEM -------------------------------
if options['method'] == 'l2':
    gfu.vec.data = a.mat.Inverse(Vh.FreeDofs(), 'umfpack') * f.vec
elif options['method'] == 'emb_trefftz':
    Tgfu = TA.Inverse(inverse='umfpack') * (PPT * (f.vec - a.mat * uf))
    gfu.vec.data = PP * Tgfu + uf
else:
    raise Exception('Unknown method only l2 and emb_trefftz possible')


l2error = sqrt(Integrate((gfu)**2 * dx, mesh))
print(f"L2 Norm of {options['method']} solution: {l2error}")


# --------------------------- VISUALISATION OUTPUT ----------------------------
# Unset mesh adaptation
mesh.deformation = None

# GUI Visualization:
Draw(deformation, mesh, 'deformation')
DrawDC(lsetp1, gfu, 0, mesh, f"uh_{options['method']}")

# VTKs
gfw = GridFunction(HDiv(mesh, order=order))
gfw.Set(w)
vtk = VTKOutput(
    mesh, coefs=[gfu, lsetp1, CF((gfw[0], gfw[1], 0))],
    names=[f"sol_{options['method']}", 'lset', 'vel'],
    filename=f"{options['method']}_convdiff_{order}_N{options['n_ref']}",
    subdivision=order, floatsize='single', legacy=False)
vtk2 = VTKOutput(
    mesh, coefs=[CF((deformation[0], deformation[1], 0))],
    names=['deformation'],
    filename=f"convdiff_deform_{order}_N{options['n_ref']}",
    subdivision=order, floatsize='single', legacy=False)
vtk2.Do()
mesh.SetDeformation(deformation)
vtk.Do()
mesh.UnsetDeformation()
Redraw()
