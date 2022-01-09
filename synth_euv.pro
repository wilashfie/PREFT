function synth_euv, tube, flux=flux, wave=wave
;  wave = [ 'xxx', '304', '171', '193', '211', '335', '94', '131' ]
;             0      1      2      3      4      5      6     7

if( not keyword_set( flux ) ) then flux = 1.0e3;  in untirs of 1.0e16 Mx

axis = tube.x
ang = 2.0/!radeg
axis[1,*] = cos(ang)*tube.x[2,*]
axis[2,*] = sin(ang)*tube.x[2,*]

rad = sqrt( flux/tube.b/!pi )

n_e = tube.epamu*(tube.rho*1.0d-16)/1.67d-24;  cm^{-3}

restore, 'aia_resp.sav'
rsp = aia_resp.a171.tresp
if( keyword_set( wave ) ) then begin
  case wave of
	1: rsp = aia_resp.a304.tresp
	2: rsp = aia_resp.a171.tresp
	3: rsp = aia_resp.a193.tresp
	4: rsp = aia_resp.a211.tresp
	5: rsp = aia_resp.a335.tresp
	6: rsp = aia_resp.a94.tresp
	7: rsp = aia_resp.a131.tresp
  endcase
endif
g = interpol( rsp, aia_resp.logte, alog10( tube.t )+6.0 )

em = n_e^2*g*1.0d8

; plot, axis[0,*], axis[1,*], xr=[-10,10]
; show_loop_image, axis, rad, em, /over, pix=0.25
; oplot, axis[0,*], axis[1,*]

xr = [ -5,5]
zr = [-5,15]
pix = 0.2

nx = ceil( (xr[1]-xr[0])/pix )
nz = ceil( (zr[1]-zr[0])/pix )

x = xr[0] + pix*findgen(nx)
z = zr[0] + pix*findgen(nz)

img = loop_image( axis, rad, em, x, z )

res = { img:img, x:x, z:z }

return, res
end
