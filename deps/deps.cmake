# set(SNOPT_SRCS
#     deps/snopt/np02lib.f
#     deps/snopt/npopt.f
#     deps/snopt/sn02lib.f
#     deps/snopt/sn03prnt.f
#     deps/snopt/sn05wrpa.f
#     deps/snopt/sn05wrpb.f
#     deps/snopt/sn05wrpc.f
#     deps/snopt/sn05wrpn.f
#     deps/snopt/sn10mach.f
#     deps/snopt/sn12ampl.f
#     deps/snopt/sn15blas.f
#     deps/snopt/sn17util.f
#     deps/snopt/sn20amat.f
#     deps/snopt/sn25bfac.f
#     deps/snopt/sn27lu.f
#     deps/snopt/sn30spec.f
#     deps/snopt/sn35mps.f
#     deps/snopt/sn37wrap.f
#     deps/snopt/sn40bfil.f
#     deps/snopt/sn50lp.f
#     deps/snopt/sn55qp.f
#     deps/snopt/sn56qncg.f
#     deps/snopt/sn57qopt.f
#     deps/snopt/sn60srch.f
#     deps/snopt/sn65rmod.f
#     deps/snopt/sn70nobj.f
#     deps/snopt/sn80ncon.f
#     deps/snopt/sn85hess.f
#     deps/snopt/sn87sopt.f
#     deps/snopt/sn90lmqn.f
#     deps/snopt/sn95fmqn.f
#     deps/snopt/snopta.f
#     deps/snopt/snoptb.f
#     deps/snopt/snoptc.f
#     deps/snopt/snoptq.f
#     deps/snopt/sq02lib.f
#     deps/snopt/sqopt.f)
#
# set(ORDERPACK_SRCS
#     deps/orderpack/ctrper.f90
#     deps/orderpack/fndnth.f90
#     deps/orderpack/givcor.f90
#     deps/orderpack/indmed.f90
#     deps/orderpack/indnth.f90
#     deps/orderpack/inspar.f90
#     deps/orderpack/inssor.f90
#     deps/orderpack/median.f90
#     deps/orderpack/mrgref.f90
#     deps/orderpack/mrgrnk.f90
#     deps/orderpack/mulcnt.f90
#     deps/orderpack/rapknr.f90
#     deps/orderpack/refpar.f90
#     deps/orderpack/refsor.f90
#     deps/orderpack/rinpar.f90
#     deps/orderpack/rnkpar.f90
#     deps/orderpack/uniinv.f90
#     deps/orderpack/unipar.f90
#     deps/orderpack/unirnk.f90
#     deps/orderpack/unista.f90
#     deps/orderpack/valmed.f90
#     deps/orderpack/valnth.f90)
#
# add_library(snopt ${SNOPT_SRCS})
# add_library(orderpack ${ORDERPACK_SRCS})
#
#
# OPTION(FoX_ENABLE_WCML OFF)
# OPTION(FoX_ENABLE_WKML OFF)
# OPTION(FoX_ENABLE_EXAMPLES OFF)
# add_subdirectory(deps/fox)
