pro set_alfven_energy_densities, tube
; old
; now just set to zero.
;
; tube made, bent, now adding new variables from preft_awave:
; take care of Elsasser variables first, so we can deal with scalar energy densites from here on.

tube.va = tube.b/sqrt(4*!pi*tube.rho_e)

zp = tube.v + ([1,1,1]#tube.va)
zm = tube.v - ([1,1,1]#tube.va)

; dot and center Elsasser arrays
zp2 = total(zp*zp, 1)
zp2c = 0.5*( shift( zp2, -1 ) + zp2 )
zm2 = total(zm*zm, 1)
zm2c = 0.5*( shift( zm2, -1 ) + zm2 )

tube.wp = tube.rho*zp2c/4.0
tube.wm = tube.rho*zm2c/4.0

return
end
