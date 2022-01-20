plt = 0 ; boolean for plotting

l_fin = 60.0
tmax = 4.0 ;MK
tmin = 0.01
dth = 120.0  ;  theta0 = 100.0
b0 = 100.0

derg = 6.2d11 ; erg/cm2

l_rx = 1.0d-8*2*!pi*derg/(b0*sin(0.25*dth/!radeg)^2)^2;  length of retracted section
d_len = 2*l_rx*sin(0.25*dth/!radeg)^2;   change due to retaction
l0 = l_fin + d_len


init_lam_rlf, tmin=1.0d6*(1.01*tmin)
preft_awave, /verb
set_uniform_cs, b0
tube = rtv_tube_equilib( l0, tmin=tmin, chr=3.0, tmax=tmax, n=200 )
bend_tube, tube, dth;, bend_frac=0.25; , /zero
set_heat_for_equilib, tube
print, 'n=', tube.n
print, 'min(dl) = ', min(tube.dl)


; ===== parameters for Drag
;tube.drag_params[0] = [ 0.2 ] ; frac x-ferred to heat
tube.drag_params[0] = [ 0.0 ] ; frac x-ferred to heat
; ==== the drag
;drag_coeff = 0.0;  in Mm^{-1}
drag_coeff = 8.0;  in Mm^{-1}
min_z = min( tube.x[2,*] );  feet
;  no drag @ feet
dz_drag_free = 8.0;  Mm
tube.drag_const = drag_coeff*( 0.5 + 0.5*tanh( 5*( tube.x[2,*] - min_z + dz_drag_free) ) )

; ====== parameters 4 Alfven energy densities
;tube.alfven_params[0:2] = [0.05, 0.1, 0.1, 0.95] ; heat, kperp, kappa adjust, reflection, p-wave
;tube.alfven_params[0:4] = [0.0, 5.0, 1.0, 0.5, 0.0] ;  for test.
tube.alfven_params[0:4] = [0.5, 9.0, 1.0, 0.5, 1.0] ;  for test.

; nte
init_nte_heat, tube, delta=5.0, alpha=4.0, ec=20.0 ; power-law, beaming parameter, cut-off in keV,   (?) #different flux tube segments


; for retraction time
va = b0/sqrt(4*!pi*tube.rho[tube.n/2])
t_rx = 0.5*l_rx/va;  time to permit reconnection

calc_len,tube
print, 'tube length = ', l0
print, '  reconn. length = ', l_rx
print, 'reconnection time = ', t_rx



tube.inv_hflf = 1.0




; now the run
dt = 0.1
nst = ceil( t_rx/dt )
tarr = [ tube ]


; plotting stuff
if (plt eq 1) then begin
window, 0, xs=600, ys=500
plot, tube.x[0,*], tube.x[2,*], /iso, yr=[-30,2], yst=1

window, 1, xs=600, ys=500
window, 2, xs=600, ys=900
window, 3, xs=600, ys=500
endif


i = 1
repeat begin
  adv_tube, tube, dt, max=500000L
  calc_len, tube
  l_cur = max( tube.l )
  tarr = [ tarr, tube ]
  print, 'Tmax = ', max(tube.t)
  if(plt eq 1) then begin
     wset, 1
     plot, tube.x[0,*], tube.x[2,*], /iso, yr=[-30,2], yst=1, tit=string(tube.time)
     ;plot_io, tube.l, tube.p
     wset, 2
     tube_plot, tube, fp=5.0
     wset, 3
     plot, tube.l, tube.wm+tube.wp;, psym = 1;, yr = [0.0,500]
     oplot, tube.l, tube.wm, color = cgcolor('blue');, psym = 1;, yr = [0.0,500]
     oplot, tube.l, tube.wp, color = cgcolor('red');, psym = 1;, yr = [0.0,500]
    ;plot, tube.l, tube.kap
    ;oplot, tube.l, tube.kap_old, color = cgcolor('red')
    ;oplot, tube.l, shift( tube.wp, -1 ) - tube.wp, psym = 2, color = cgcolor('red');, yr = [0.0,500]
    ;oplot, tube.l, tube.wm, color = cgcolor('red')
  endif
  save, file='~/data/drag/awave/turb2.sav', tarr
endrep until( (l_cur lt l_fin) )


stop


; now straighten and continue run
straighten_tube, tube, cf_smooth=0.01
dt = 0.2
nst = 300

;window, 1, xs=600, ys=900
;window, 2, xs=600, ys=500
for i=1, nst do begin
  adv_tube, tube, dt, max=500000L
  tarr = [ tarr, tube ]
;  wset, 1
;  tube_plot, tube, fp=7.0, vfp=0.25
;  print, tube.time, max( tube.t )
;  save, file='~/data/drag/apr_flare_cur.sav', tarr
endfor


;  prepare energy array
erg = [ tube_erg( tarr[0] ) ]
for j=1, nst do erg = [ erg, tube_erg( tarr[j] ) ]

end
