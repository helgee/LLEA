program testpropgators

use types, only: dp
use assertions
use constants
use epochs
use exceptions
use math, only: eps
use propagators
use states
use trajectories

implicit none

type(exception) :: err
real(dp) :: mu
real(dp) :: dt
real(dp), dimension(6) :: rv0exp
real(dp), dimension(6) :: rv1exp
real(dp), dimension(6) :: rv0
real(dp), dimension(6) :: rv1
type(kepler) :: kep
type(epoch) :: ep
type(state) :: s0
type(state) :: s1
type(epochdelta) :: epd
type(trajectory) :: tra

call init_constants

! Source: Vallado, Fundamentals of Astrodynamics and Applications, 4th edition, p. 94-95
mu = 3.986004418e5_dp
rv0exp = [1131.340_dp, -2282.343_dp, 6672.423_dp, -5.64305_dp, 4.30333_dp, 2.42879_dp]
dt = 40._dp * 60._dp
rv1exp = [-4219.752737795695_dp, 4363.029177180833_dp, -3958.7666166029762_dp, &
    3.689866025052519_dp, -1.9167347770873095_dp, -6.112511100000719_dp]
rv1 = solve_kepler(mu, rv0exp, dt)
call assert_almost_equal(rv1, rv1exp, __LINE__)
rv0 = solve_kepler(mu, rv1, -dt)
call assert_almost_equal(rv0, rv0exp, __LINE__)
rv1 = solve_kepler(mu, rv0exp, eps)
call assert_almost_equal(rv1, rv0exp, __LINE__)
rv1 = solve_kepler(mu, rv0exp, dt, iterations=1, err=err)
call assert_raises("Kepler solver did not converge.", err, __LINE__)
call reset(err)

ep = epoch()
s0 = state(ep, rv0exp)
s1 = getstate(s0, dt, kep)
call assert_almost_equal(s1%rv, rv1exp, __LINE__)
epd = epochdelta(seconds=dt)
s1 = getstate(s0, epd, kep)
call assert_almost_equal(s1%rv, rv1exp, __LINE__)
kep = kepler(iterations=1)
s1 = getstate(s0, epd, kep, err=err)
call assert_raises("Kepler solver did not converge.", err, __LINE__)
call reset(err)

kep = kepler()
tra = gettrajectory(s0, epd, kep)
call assert_almost_equal(tra%final_state%rv, rv1exp, __LINE__)

end program testpropgators