;***************************************************************
   function imake_scalar_accum,nstart,nstop,nstep,nfreq,magnetic,iswap=iswap,flt=flt
;---------------------------------------------------------------
; Makes accumulation of flux files
; Nrecs = ndumps per file
;
; Scalars are :
;       0 = n_step
;       1 = total_time
;       2 = u_min
;       3 = u_max,
;       4 = u_avg
;       5 = v_min
;       6 = v_max
;       7 = v_avg
;       8 = w_min
;       9 = w_max
;      10 = w_avg
;      11 = kinetic_energy
;      12 = max_kinetic_energy
;      13 = bx_min
;      14 = bx_max
;      15 = b_avg(1)
;      16 = by_min
;      17 = by_max
;      18 = b_avg(2)
;      19 = bz_min
;      20 = bz_max
;      21 = b_avg(3)
;      22 = Magnetic_Energy
;      23 = Max_MagneticPressure
;---------------------------------------------------------------

files=file_search('????????')

if n_elements(nstart) eq 0 then nstart = 0
if n_elements(nstop) eq 0 then begin
  s = size(files)
  len = s(1)
  nstop = Long(files(len-1))
endif
if n_elements(nstep) eq 0 then nstep = Long(files(0))
if n_elements(nfreq) eq 0 then nfreq = 50
if n_elements(magnetic) eq 0 then magnetic = 1

;print,nstart,nstop,nstep,nfreq,magnetic

if n_elements(flt) eq 0 then flt = 0

nf = (long(nstop) - nstart) / nstep
nrec = (long(nstop) - nstart) / nfreq
rpf = long(nstep) / nfreq

if(magnetic NE 0) then begin
  ns = 24
endif else begin
  ns = 13
endelse

if flt eq 0 then begin
  a = dblarr(ns,nrec)
endif else begin
  a = fltarr(ns,nrec)
endelse

cur = 0

if flt eq 0 then begin
  in = assoc(1, dblarr(ns, rpf)) 
endif else begin
  in = assoc(1, fltarr(ns, rpf)) 
endelse

for i=long(nstep),long(nstop),long(nstep) do begin
  file = string(i,format='(I8.8)')
  openr,1,file
  temp = in(0)
  close,1
  a(*,cur:cur+rpf-1) = temp(*,*)
cur = cur+rpf
endfor

return,a

end
