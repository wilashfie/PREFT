; report the characteristics of this run based on its initial state

pro bend_smry, tarr, unit=unit

tube = tarr[0]
calc_len, tube

n = tube.n
nh = n/2

b = tube.b[nh]
tmax = tube.t[nh]
ne_min = tube.rho[nh]*tube.epamu/1.67d-8
tmin = tube.t[3]
ii = where( tube.t gt 1.1*tmin )
l0 = min( tube.l[ii] )
l1 = max( tube.l[ii] )

tv1 = tube.x[*,20] - tube.x[*,10]
tv1m = sqrt( total( tv1^2 ) )
tv2 = tube.x[*,n-10] - tube.x[*,n-20]
tv2m = sqrt( total( tv2^2 ) )
cang = total( tv1*tv2 )/tv1m/tv2m
ang = !radeg*acos( cang )

i0 = min( ii )
cor_col = 0.5*total( tube.rho[ii]*tube.dl[ii] )*tube.epamu/1.67d-16
chr_col = total( tube.rho[0:i0]*tube.dl[0:10] )*tube.epamu/1.67d-16

nop = 11
op_arr = replicate( '   ', nop )
op_arr[0] = ' n = ' + string( n )
op_arr[1] = ' ang = ' + string( ang )
op_arr[2] = ' b = ' + string( b )
op_arr[3] =  ' apex: T = ' + string( tmax ) +  ',   n_e = ' + string( ne_min )
op_arr[4] = ' length = ' + string( l1-l0 )
op_arr[5] = ' coronal e column = ' + string( cor_col ) + ' cm^(-2)'
op_arr[6] = ' chromosphere:'
op_arr[7] = '  T =' + string( tmin )
op_arr[8] = '  l =' + string( l0 ) + ', and ' + string( max( tube.l)-l1 )
op_arr[9] = '  electron column = ' + string( chr_col ) + ' cm^(-2)'

; time until straightened
x0 = reform( tarr.x[0,0] )
d = shift( x0, -1 ) - x0
ii = where( abs(d) gt 0.1, nii )
if( nii lt 1 ) then begin
  op_arr[10] = 'never straightened'
endif else begin
  op_arr[10] ='straightened @ t=' + string( tarr[ii[0]].time )
endelse

if( keyword_set( unit ) ) then for i=0, nop-1 do printf, unit, op_arr[i] $
  else for i=0, nop-1 do print, op_arr[i]

return
end
