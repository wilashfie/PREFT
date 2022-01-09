lt0 = [  4.000, 4.896, 5.419, 5.563, 6.183, 6.563, 6.978, 7.467 ]
lt1 = [  4.896, 5.419, 5.563, 6.183, 6.563, 6.978, 7.467, 9.000 ]
lf0 = [  -29.411, -21.927, -10.565, -22.849,  -8.679, -23.867, -13.248, -25.105 ]
p = [    1.659,   0.131,  -1.966,   0.242,  -2.050,   0.264,  -1.257,   0.331 ]

print, '\Lambda(T) = \left\{ \begin{array}{lrcl}'

nr =  n_elements( lt0 )
FOR i=0, nr-1 DO BEGIN
  print, '10^{', string( lf0[i], form='(f7.2)' ), '} ~T^{'
  print, string( p[i], form='(f6.2)' ), '}~~,~~&'
  if( i gt 0 ) then $
    print, '10^{', string( lt0[i], form='(f5.2)' ), '} <'
  print, '& T &'
  if( i lt nr-1 ) then $
  print, '< 10^{', string( lt1[i], form='(f5.2)' ), '}\\'
endfor

print, '\end{array}\right.'

end