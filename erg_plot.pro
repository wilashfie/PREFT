; assuming everything is in tarr

n = n_elements( tarr )
erg = tube_erg( tarr[0] )
for i=1, n-1 do erg = [ erg, tube_erg( tarr[i] ) ]

dx = tarr[0].x[*,tarr[0].n-1] -  tarr[0].x[*,0]
me0 = sqrt( total( dx^2 ) )*tarr[0].b[0]*0.25/!pi

plot, tarr.time, erg.mag-me0, pos=[0.1, 0.1, 0.9, 0.49], /nodata, $
  ytit='10!u27!n erg'
oplot, tarr.time, erg.mag-me0, color=2
oplot, tarr.time, erg.tot-me0, lines=2
oplot, tarr.time, erg.kin, color=4
oplot, tarr.time, erg.kin_par, color=3
ylab=0.44
xyouts, 0.2, ylab, 'free mag.', /norm, color=2
xyouts, 0.35, ylab, 'kinetic', /norm, color=4
xyouts, 0.5, ylab, 'kinetic ||', /norm, color=3

; blow-up
plot, tarr.time, erg.therm-erg[0].therm, pos=[0.1, 0.51, 0.9, 0.9], /nodata, $
  /noerase, xtickf='no_tick_label',  ytit='10!u27!n erg', $
  tit='Energies for 10!u19!n Mx tube'
oplot, tarr.time, erg.kin_par, color=3
oplot, tarr.time, tarr.net_erg_loss/tarr[0].b, color=6
oplot, tarr.time, erg.therm-erg[0].therm, color=8

ylab=0.84
xyouts, 0.2, ylab, 'kinetic ||', /norm, color=3
xyouts, 0.35, ylab, 'thermal', /norm, color=8
xyouts, 0.5, ylab, 'radiated', /norm, color=6


end