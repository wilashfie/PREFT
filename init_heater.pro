; a set of programs which will update the heating profile of a tube
; over time.

pro apply_heater, tube

common heater_params, h0, hprof, dur_h, hmax

if( tube.time gt 1.1*dur_h ) then return

amp = ( 1.0 - 2*abs( tube.time - 0.5*dur_h )/dur_h ) > 0.0
tube.heat = h0 + amp*hmax*hprof

return
end

; ---------------------------------------------

pro init_heater, tube, etot, dur=dur, len=len, l0=l0

common heater_params, h0, hprof, dur_h, hmax

h0 = tube.heat;  save the initial heat to be augmented
calc_len, tube

if( not keyword_set( l0 ) ) then l0 = 0.5*max( tube.l );  midpoint
if( not keyword_set( len ) ) then len = 0.5*max( tube.l );  heat center 50%
if( not keyword_set( dur ) ) then dur = 10.0;  10 seconds of heating
dur_h = dur

hprof = ( 1.0 - 2*abs( tube.l - l0 )/len ) > 0.0;  normalized tent
hprof = hprof*tube.dl;   will store heat in each cell

htot = total( hprof )
hprof = hprof/htot;  normalized: integrates to 1
hmax = 2*etot/dur;  for tent time profile

return
end



