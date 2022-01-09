; a set of programs which will update the heating profile of a tube
; over time.

pro apply_nte_source, tube

common nte_params, h0, hprof, dur_h, hmax

if( tube.time gt 1.1*dur_h ) then return

amp = ( 1.0 - 2*abs( tube.time - 0.5*dur_h )/dur_h ) > 0.0
tube.heat = h0 + amp*hmax*hprof

return
end

; =====================================

function psi_func_calc, x, alpha=alpha

if( n_elements( alpha ) lt 1 ) then alpha = 5.0

n = n_elements( x )

dx = shift( x, -1 ) - x
xm = 0.5*( shift( x, -1 ) + x )
xm[n-1] = x[n-1]

ii = where( xm lt 1.0, nii )

n = max(ii)

f = ( 1.0 - xm ) > 0.0
f = xm^(-alpha) * f^(-0.666667)

; the x ~ 1 limit
psi = ( 1.0 - x ) > 0.0
psi = 3.0*psi^(0.333333)

for i=6, n do psi[n-i] = psi[n-i+1] + f[n-i]*dx[n-i]

psi = psi*x^(alpha-1.0)


return, psi
end

; =====================================

function w_func_calc, xi_max, alpha=alpha, delta=delta

ngrid = 300
x = xi_max*( findgen( ngrid ) + 1.0 )/float(ngrid);  create grid

if( n_elements( alpha ) lt 1 ) then alpha = 5.0
if( n_elements( delta ) lt 1 ) then delta = 5.0

xf = x > 0.5
q = xf^(-0.5*delta)*beta(0.5*delta,0.33333)

ii = where( x lt 1.0 )
q[ii] = x[ii]^(-0.5*delta)*ibeta( 0.5*delta, 0.33333, x[ii] )*beta(0.5*delta,0.33333)

fctr = 1.0
if( alpha le 1.0 ) then begin
   q[ii] = q[ii] + x[ii]^(alpha-1.0)* $
       ibeta( 0.33333, 1.0-alpha, 1.0 - x[ii] )*beta(0.33333,1.0-alpha)
  fctr = 2*alpha/(2*alpha+delta-2.0)
endif else if( alpha lt 10.0 ) then begin
  psi = psi_func_calc( x, alpha=alpha )
  q[ii] = q[ii] + psi[ii]
  fctr = 2*alpha/(2*alpha+delta-2.0)
endif

w = q*fctr*(delta-2.0)/6.0

; the auxiliary integrals
y0 = fltarr( ngrid )
y1 = fltarr( ngrid )

for i=1, ngrid-1
  dxi = x[i] - x[i-1]
  y0[i] = y0[i-1] + 0.5*( w[i] + w[i-1] )*dxi
  y1[i] = y0[i-1] + 0.5*( w[i]*x[i] + w[i-1]*x[i-1] )*dxi
endfor

res = { xi:x, w:w, y0:y0, y1:y1 }

return, res
end

; =====================================
;  tube.drag_heat currently contains the energy lost to drag force in [ 1.0e8 erg/sec/cm^2 ].
;  break this into a direct contribution from heating (fraction = tube.drag_params[1])
;  and that contributed to NTe's.  update tube.drag_heat and tube.nte_heat

pro calc_nte_heat, tube

common nte_params, kmat, ns, n_co, iseg, dxi, xi_mid, eta, q_wt


f_dir = tube.drag_params[1];  fraction of energy lost which is deposited directly
f_nte = 1.0 - f_dir

tube.nte_heat = f_nte*tube.drag_heat;  fraction converted to NTE
for j=0, ns-1 do eta[j] = total( tube.nte_heat[ (iseg[0,j]) : (iseg[1,j]) ] )/dxi[j]
qdot = 0.5*n_co*( kmat # eta )
;   now update the heating profile
for j=0, ns-1 do begin
  i0 = iseg[0,j]
  i1 = iseg[1,j]
  tube.nte_heat[i0:i1] = qdot[j]*q_wt[i0:i1]
endfor

tube.drag_heat = f_dir*tube.drag_heat;   allow for some fraction of direct deposition

return
end

; =====================================

pro init_nte_heat, tube, nseg=nseg, delta=delta, alpha=alpha, ec=ec

common nte_params, kmat, ns, n_co, iseg, dxi, xi_mid, eta, q_wt

calc_len, tube

if( not keyword_set( delta ) ) then delta = 4.0; the power law
if( not keyword_set( alpha ) ) then alpha = 4.0; beaming parameter (1 = isotropic)
if( not keyword_set( ec ) ) then ec = 10.0;   cut-off in keV
if( not keyword_set( nseg ) ) then nseg = 40;  number of different flux tube seg
ns = nseg;  to store in common block

; compute normalized column
rho2ne = tube.epamu/1.0d16/1.67d-24
n_e = tube.rho*rho2ne
dn = tube.dm*tube.b*rho2ne*1.0d8
ncol = dblarr( tube.n )
for i=1, tube.n-1 do ncol[i] = ncol[i-1] + dn[i-1]
n_co = 1.0d17*(ec)^2;  stopping column for cut-off
xi = ncol/n_co;  normalized column

; divide into nseg sections
dxi0 = max( xi )/float( nseg )
iseg = lonarr( 2, nseg )
iseg[0,0] = 0;              first segment starts at i=0
iseg[1,nseg-1] = tube.n-1;   last segment 
for j = 0, nseg-2 do begin
  i0 = iseg[0,j]
  mm = min( abs( xi - (j+1)*dxi0 ), i1 )
  iseg[1,j] = i1 > (i0+1)
  iseg[0,j+1] = i1+1
endfor

;  compute positions and sizes of segments
dxi = fltarr( nseg )
xi_mid = fltarr( nseg )
q_wt = fltarr( tube.n );     returns fractions of dot_j to tube cells
for j = 0, nseg-1 do begin
  i0 = iseg[0,j]
  i1 = iseg[1,j]
  xi_mid[j] = 0.5*( xi[i1] + xi[i0] );   central position of segment j (normalized column)
  dxi[j] = xi[i1] - xi[i0];              total size of segment j (normalized column)
  q_wt[i0:i1] = dn[i0:i1]/n_co/dxi[j];   fractional contribution of each cell in segment j
endfor

; compute W(xi) and components of diagonal D(xi)
dxi_max = 1.05*max( dxi ) > 1.0
res = w_func_calc( dxi_max, alpha=alpha, delta=delta )

; ===== the matrix
kmat = fltarr( nseg, nseg )

; diagonal elements
for i=0, nseg-1 do begin
  y0 = interpol( res.y0, res.dxi, dxi[i] )
  y1 = interpol( res.y1, res.dxi, dxi[i] )
  kmat[i,i] = 2.0*( dxi[i]*y0[0] - y1[0] );  kmat[i,i]= D(dxi)*( dxi^2 )
endfor

; off-diagonals
fctr = alpha*(delta-2.0)*beta(0.5*delta,0.33333)/3.0/(2*alpha - delta - 2.0 ) )
for i=0, nseg-2 do begin
  for j=i+1, nseg-1
    dxi_ij = xi[j] - xi[i];  should always be positive (!!)
    if( dxi_ij lt 1.0 ) then kmat[i,j] = interpol( res.w, res.dxi, dxi_ij ) $
      else kmat[i,j] = fctr*( dxi_ij )^(-0.5*delta )
    kmat[i,j] = kmat[i,j]*dxi[i]*dxi[j]
    kmat[j,i] = kmat[i,j];  lower triangle
  endfor
endfor

return
end



