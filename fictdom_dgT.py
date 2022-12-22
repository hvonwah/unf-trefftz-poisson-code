# ------------------------------ LOAD LIBRARIES -------------------------------
from netgen.geom2d import SplineGeometry
from netgen.csg import CSGeometry, OrthoBrick, Pnt
from ngsolve.meshes import *
from ngsolve import *
from ngsolve import dx as ngsdx
from ngsolve.internal import *
from xfem import *
from xfem.lsetcurv import *
from ngstrefftz import *
import time

ngsglobals.msg_level = 2


def SolveFdT(example=1, method='emb_trefftz', stabil='global', maxh=0.5, order=3, quad_mesh=False, struc_mesh=False, deformation=True):
    '''
    Solve a Laplace problem on the domain on a ring shaped domain using
    an unfitted (Trefftz) discontinuous Galerin method.s

    Parameters
    ----------
    example : int
        The numerical example to run. Possible are 1, 2, 3. Examples 1
        and 3 are two dimensional and example 2 is three dimensional.
    method : str
        The method to use. Can be 'emb_trefftz' for the embedded Treftz
        method, 'trefftz' for the strong Trefftz method, 'trefftz_mixed'
        for the mixed Trefftz method or 'l2' for the standard dG method.
    stabil : str
        Either 'global' or 'patch'. Whether to use global or patch-wise
        ghost-penalties. The patch-wise operator implements the "weak"
        aggregated methods.
    maxh : float
        Mesh size of the domain.
    order: int
        Polynomial order of the finite element space.
    quad_mesh : bool
        Construct a quadrilateral mesh.
    struc_mesh : bool
        Construct a structured mesh of the domain.
    deformation : bool
        Deform the mesh to realise high-order geometry approximation.

    Returns
    -------
    tuple(l2 err, ndof)
        Tuple with the L2 error and the numbers of degrees of freedom
        on which the problem is solved.
    '''

    if method not in ['l2', 'emb_trefftz', 'trefftz', 'trefftz_mixed']:
        raise Exception('Unknown method passed to SolveFdT')

    # Stabilization parameter for ghost-penalty
    gamma_stab = 0.01
    # Stabilization parameter for Nitsche
    lambda_nitsche = 10 * order * order
    # Stabilization parameter for interior penalty
    lambda_dg = 10 * order * order

    # -------------------------- GEOMETRY AND MESH ----------------------------
    if example in [1, 3] and struc_mesh is False:
        square = SplineGeometry()
        square.AddRectangle((-1, -1), (1, 1), bc=1)

        ngmesh = square.GenerateMesh(maxh=maxh, quad_dominated=quad_mesh)
        mesh = Mesh(ngmesh)
    elif example in [1, 3] and struc_mesh is True:
        print('Using strctured 2d mesh')
        space_refs = (int)(log(maxh) / log(0.5))
        print('ref exponent: ', space_refs)
        mesh = MakeStructured2DMesh(quads=quad_mesh, nx=2**(space_refs + 1),
                                    ny=2**(space_refs + 1),
                                    mapping=lambda x, y: (2 * x - 1, 2 * y - 1))
    elif example == 2 and struc_mesh is False:
        cube = CSGeometry()
        cube.Add(OrthoBrick(Pnt(-1, -1, -1), Pnt(1, 1, 1)))
        mesh = Mesh(cube.GenerateMesh(maxh=maxh))
    elif example == 2 and struc_mesh is True:
        if quad_mesh:
            print('WARNING: CutDG with hexes is experimental')
        print('Using structured 3d mesh')
        space_refs = (int)(log(maxh) / log(0.5))
        print('ref exponent: ', space_refs)
        mesh = MakeStructured3DMesh(hexes=quad_mesh, nx=2**(space_refs + 1),
                                    mapping=lambda x, y, z: (2 * x - 1,
                                                             2 * y - 1,
                                                             2 * z - 1))

    # ---------------------------- EXACT SOLUTION -----------------------------
    if example in [1, 3]:
        r2 = 3 / 4  # outer radius
        r1 = 1 / 4  # inner radius
        rc = (r1 + r2) / 2.0
        rr = (r2 - r1) / 2.0
        r = sqrt(x**2 + y**2)
        levelset = IfPos(r - rc, r - rc - rr, rc - r - rr)
    elif example == 2:
        levelset = sqrt(x**2 + y**2 + z**2) - 0.5
        levelset += 0.5 / 3.5 * cos(5 * atan2(y, x)) * cos(pi * z)
    else:
        raise ValueError('Unknow Example requested')

    if example == 1:
        exact = (exp(x) * sin(y)).Compile()
        coef_f = 0
        coef_g = exact
    elif example == 2:
        exact = (exp(sqrt(2) * x) * sin(y) * cos(z)).Compile()
        coef_f = 0
        coef_g = exact
    elif example == 3:
        exact = (20 * (r2 - r) * (r - r1)).Compile()
        coef_f = - (exact.Diff(x).Diff(x) + exact.Diff(y).Diff(y)).Compile()
        coef_g = 0
        if method == 'trefftz':
            raise Exception('Trefftz only works for homogeneous data')
    # -------------------- LEVEL SET AND H-O APPROXIMATION --------------------
    lsetmeshadap = LevelSetMeshAdaptation(mesh, order=order, threshold=0.1,
                                          discontinuous_qn=True,
                                          heapsize=int(1e8))
    if deformation is True:
        print('Deforming mesh')
        deformation = lsetmeshadap.CalcDeformation(levelset)
        mesh.SetDeformation(deformation)
    else:
        InterpolateToP1(levelset, lsetmeshadap.lset_p1,
                        eps_perturbation=lsetmeshadap.eps_perturbation)
    lsetp1 = lsetmeshadap.lset_p1

    # ----------------------- ELEMENT AND FACET MARKERS -----------------------
    ci = CutInfo(mesh, lsetp1)
    els_hasneg = ci.GetElementsOfType(HASNEG)
    els_if = ci.GetElementsOfType(IF)
    # Facets with parts inside the domain
    facets_dg = GetFacetsWithNeighborTypes(mesh, a=els_hasneg, b=els_hasneg)
    # Facets for ghost penalty stabilization
    if stabil == 'patch':
        els_neg = ci.GetElementsOfType(NEG)
        EA = ElementAggregation(mesh)
        EA.Update(els_neg, els_if)
        facets_gp = EA.GetInnerPatchFacets()
    else:  # if stabil == 'global':
        facets_gp = GetFacetsWithNeighborTypes(mesh, a=els_hasneg, b=els_if)

    # ------------------------- FINITE ELEMENT SPACE --------------------------
    if method == 'trefftz_mixed':
        Vhbase_1 = L2(mesh, order=order, dirichlet=[], dgjumps=True)
        Vh_1 = Restrict(Vhbase_1, els_hasneg)

        Vhbase_om2 = L2(mesh, order=order - 2, dirichlet=[], dgjumps=True)
        Vh_om2 = Restrict(Vhbase_om2, els_hasneg)

        Vh = Vh_1 * Vh_om2
        gfU = GridFunction(Vh)
        gfu, gfp = gfU.components
        (u, p), (v, q) = Vh.TnT()
    elif method == 'trefftz':
        Vhbase = trefftzfespace(mesh, order=order, eq='laplace', dgjumps=True)
        Vh = Restrict(Vhbase, els_hasneg)
        u, v = Vh.TnT()
        gfu = GridFunction(Vh)
        gfU = gfu
    elif method == 'emb_trefftz':
        Vhbase = monomialfespace(mesh, order=order, dirichlet=[], dgjumps=True)
        Vh = Restrict(Vhbase, els_hasneg)
        # Vhbase_om2 = L2(mesh, order=order - 2, dirichlet=[], dgjumps=True)
        Vhbase_om2 = monomialfespace(mesh, order=order - 2, dirichlet=[], dgjumps=True)
        Vh_om2 = Restrict(Vhbase_om2, els_hasneg)
        u, v = Vh.TnT()
        gfu = GridFunction(Vh)
        gfU = gfu
    else:
        Vhbase = monomialfespace(mesh, order=order, dirichlet=[], dgjumps=True)
        Vh = Restrict(Vhbase, els_hasneg)
        u, v = Vh.TnT()
        gfu = GridFunction(Vh)
        gfU = gfu

    # ------------------------------ INTEGRATORS ------------------------------
    h = specialcf.mesh_size
    nh = 1.0 / Norm(grad(lsetp1)) * grad(lsetp1)
    nF = specialcf.normal(mesh.dim)
    flux_u = -0.5 * (grad(u) + grad(u.Other())) * nF
    flux_v = -0.5 * (grad(v) + grad(v.Other())) * nF
    jump_u = u - u.Other()
    jump_v = v - v.Other()

    dX = ngsdx(definedonelements=els_hasneg)
    dx = dCut(lsetp1, NEG, definedonelements=els_hasneg)
    dk = dCut(lsetp1, NEG, skeleton=True, definedonelements=facets_dg)
    ds = dCut(lsetp1, IF, definedonelements=els_if)
    dw = dFacetPatch(definedonelements=facets_gp)

    def Lap(u):
        return sum(Trace(u.Operator('hesse')))

    a = RestrictedBilinearForm(Vh, element_restrcition=els_hasneg,
                               facet_restriction=facets_dg,
                               check_unused=False)
    a += (grad(u) * grad(v)).Compile() * dx
    a += (lambda_dg / h * jump_u * jump_v + flux_u * jump_v + flux_v * jump_u).Compile() * dk
    a += (-grad(u) * nh * v - grad(v) * nh * u + lambda_nitsche / h * u * v).Compile() * ds
    if stabil in ['patch','global']:
        a += (gamma_stab / h**2 * (u - u.Other()) * (v - v.Other())).Compile() * dw
    if method == 'trefftz_mixed':
        a += (-Lap(u) * q - Lap(v) * p).Compile() * dX
        #a += p * q * dX

    f = LinearForm(Vh)
    f += (coef_f * v).Compile() * dx
    f += (- grad(v) * nh * coef_g + lambda_nitsche / h * coef_g * v).Compile() * ds
    if method == 'trefftz_mixed':
        f += (coef_f * q).Compile() * dX

    # -------------------------- ASSEMBLE AND SOLVE ---------------------------
    times = {}

    tic = time.perf_counter()
    f.Assemble()
    a.Assemble()
    times['assemble'] = time.perf_counter() - tic

    ncutels = sum(els_hasneg)
    print(f'n cut els: {ncutels}')

    if method == 'emb_trefftz':
        tic = time.perf_counter()
        eps = 10**-7
        vt = Vh_om2.TestFunction()
        op = Lap(u) * vt * dX \
             # (gamma_stab / h**2 * (u - u.Other()) * (vt - vt.Other())).Compile() * dw
        lop = -coef_f * vt * dX
        PP, uf = TrefftzEmbedding(bf=op, fes=Vh, lf=lop, eps=eps, test_fes=Vh_om2)
        print(f'PP.shape: {PP.shape} per element: {PP.shape[0] / ncutels} -> {PP.shape[1] / ncutels}')
        PPT = PP.CreateTranspose()
        TA = PPT @ a.mat @ PP
        times['embedding'] = time.perf_counter() - tic

        tic = time.perf_counter()
        Tgfu = TA.Inverse(inverse='sparsecholesky') * (PPT * (f.vec - a.mat * uf))
        gfu.vec.data = PP * Tgfu + uf
        times['solve'] = time.perf_counter() - tic

        nd = PP.shape[1]
    elif method == 'trefftz_mixed':
        tic = time.perf_counter()
        gfU.vec.data = a.mat.Inverse(Vh.FreeDofs()) * f.vec
        times['solve'] = time.perf_counter() - tic
        times['embedding'] = 0.0

        nd = Vh_1.ndof - Vh_om2.ndof
    elif method in ['l2', 'trefftz']:
        nd = Vh.ndof
        tic = time.perf_counter()
        gfu.vec.data = a.mat.Inverse(Vh.FreeDofs(), 'sparsecholesky') * f.vec
        times['solve'] = time.perf_counter() - tic
        times['embedding'] = 0.0

    # ----------------------------- COMPUTE ERROR -----------------------------
    l2error = sqrt(Integrate((gfu - exact)**2 * dx, mesh))
    print(f'Method {method} + {stabil}, L2 Error: {l2error} ndof/el: {nd / ncutels}')
    # import netgen.gui
    # Draw(exact,mesh,'sol', min=0, max=0.25)
    # DrawDC(levelset, -1.0, 2.0, mesh, 'x')
    # DrawDC(levelset, gfu, 0, mesh, 'u', min=0, max=0.25)
    # DrawDC(levelset, exact, 0, mesh, 'exact', min=0, max=0.25)
    # input()
    return (l2error, nd, times)


if __name__ == '__main__':
    SolveFdT(method='l2', maxh=0.5)
    SolveFdT(method='l2', maxh=0.25)
    SolveFdT(method='l2', maxh=0.125)
    SolveFdT(method='l2', maxh=0.5, stabil='none')
    SolveFdT(method='l2', maxh=0.25, stabil='none')
    SolveFdT(method='l2', maxh=0.125, stabil='none')
    SolveFdT(method='trefftz', maxh=0.5)
    SolveFdT(method='trefftz', maxh=0.25)
    SolveFdT(method='trefftz', maxh=0.125)
    SolveFdT(method='trefftz_mixed', maxh=0.5)
    SolveFdT(method='trefftz_mixed', maxh=0.25)
    SolveFdT(method='trefftz_mixed', maxh=0.125)
    SolveFdT(method='emb_trefftz', maxh=0.5)
    SolveFdT(method='emb_trefftz', maxh=0.25)
    SolveFdT(method='emb_trefftz', maxh=0.125)
