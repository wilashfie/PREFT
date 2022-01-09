; flare params
l0_set = 110.0
l_fin = 85.0
tmax = 4.0
theta0 = 100.0
b0 = 100.0


;  create flux tube
tmin = 0.01
init_lam_rlf, tmin=1.0d6*(tmin*1.01)
preft_awave, /verb
set_uniform_cs, b0
tube = rtv_tube_equilib( l0_set, tmin=tmin, tmax=tmax, chr=3.0, n=200 )
bend_tube, tube, theta0, bend_frac=0.5
set_heat_for_equilib, tube
set_alfven_energy_densities, tube
print, 'n=', tube.n
print, 'min(dl) = ', min(tube.dl)

calc_len,tube
l0 = max(tube.l)

; ===== parameters for Drag
tube.drag_params[0] = [ 0.0 ] ; frac x-ferred to heat
; ==== the drag
drag_coeff = 0.5;  in Mm^{-1}
min_z = min( tube.x[2,*] );  feet
;  no drag @ feet
dz_drag_free = 8.0;  Mm
tube.drag_const = drag_coeff*( 0.5 + 0.5*tanh( 5*( tube.x[2,*] - min_z + dz_drag_free) ) )

; ====== parameters 4 Alfven energy densities
;tube.alfven_params[0:2] = [0.05, 0.1, 0.1, 0.95] ; heat, kperp, kappa adjust, reflection
tube.alfven_params[0:2] = [0.00, 0.0, 0.0, 0.0] ;  for test.
print, 'L0 = ', l0



tube.inv_hflf = 1.0



window, 0, xs=600, ys=500

plot, tube.x[0,*], tube.x[2,*], /iso, yr=[-30,2], yst=1

; now the run

window, 1, xs=600, ys=500

dt = 0.1
ttot = 5.0
nst = ceil( ttot/dt )
tarr = [ tube ]


;nst = 4 ; 4 steps, fingers crossed

window, 2, xs=600, ys=900
window, 3, xs=600, ys=500

for i=1, nst do begin
  adv_tube, tube, dt, max=500000L
  tarr = [ tarr, tube ]
  wset, 1
  plot, tube.x[0,*], tube.x[2,*], /iso, yr=[-30,2], yst=1, tit=string(tube.time)
  wset, 2
  tube_plot, tube, fp=5.0
  wset, 3
  ;plot, tube.x[0,*], tube.v[2,*]
  plot, tube.l, tube.wp
  wset, 0
  save, file='~/data/drag/awave/test5.sav', tarr
endfor

stop

;  prepare energy array
erg = [ tube_erg( tarr[0] ) ]
for j=1, nst do erg = [ erg, tube_erg( tarr[j] ) ]

end
