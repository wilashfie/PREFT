pro adv_tube_w_heat, tube, dtime, dtmax=dtmax, maxstep=maxstep, safe=safe

if( not keyword_set( dtmax ) ) then dtmax=0.01
if( not keyword_set( maxstep ) ) then maxstep=10000L
if( not keyword_set( safe ) ) then safe=0.1

dtmin=1.0d-9
t=dtmin
cnt = 0L 

; take small step to prepare for clf_step
one_step, tube, dtmin, isoth=isotherm

repeat begin
  dtc = cfl_step( tube, [1,1,1,1,0] )
  dtf = dtime - t
  dt = min( [ dtmax, dtc*safe, dtf ] ) > dtmin
  apply_heater, tube
  one_step, tube, dt, isoth=isotherm
  t = t + dt
  cnt = cnt + 1
endrep until( ( t ge dtime ) or ( cnt gt maxstep ) )

return
END
